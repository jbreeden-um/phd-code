
using DifferentialEquations
using Plots
using Printf
using LinearAlgebra

try
    Controllers
catch
    include("Dynamics.jl"); include("Constraints.jl"); include("Controllers.jl"); include("Disturbances.jl");
end

# The reader is welcome to test this code with different parameters, settings, tolerances, etc.
# The code below describes what was done to generate the plots in the paper.
# There are also a few commented out sections of code throughout this project, which should be treated as tools
# that may improve execution when settings are changed.

"""
The following function is used for the integration of only the landing simulation.
For the docking simulation equivalent, see the Controllers.jl file
"""
function DomainLimiter(x, p, t)
    return Constraints.last_H > 0 || isnan(Controllers.last_u[1]);
end

function quit_early(x, t, integrator)
    if isfile("quit.txt");
        rm("quit.txt");
        return true;
    else
        return false;
    end
end

function save_ode_data_landing(x,t,integrator)
    return [x; Controllers.last_H; Controllers.last_u[1:3]; Controllers.last_sigma; Constraints.last_t; Controllers.last_residual];
end

function save_ode_data_docking(x,t,integrator)
    return [x; Controllers.last_H_left; Controllers.last_H_right; Controllers.last_H; Controllers.last_u[1:2]; 0; Controllers.last_sigma; Constraints.last_t; Controllers.last_residual];
end

# case = :Landing;
case = :Docking;

if case==:Landing
    Controllers.set_f(Constraints.f_land);
    Controllers.set_g(Constraints.g_land);

    r0 = [-1500e3; -100e3; 50e3];
    v0 = [100.0; 0.0; 0.0];
    x0 = [r0; v0];
    u_max = 0.5;
    w_umax = 0.025;
    w_xmax = 0.01;
    gamma1 = 0.1;
    gamma2 = 1.5;
    l_h = 1;
    delta = gamma2^2/2;
    if delta < 2*(l_h*w_xmax + gamma1/2)^2
        println("Warning: Landing cannot be achieved for these tolerances");
    end
    Constraint_func = (t,x)->Constraints.H_variable(t,x,delta=delta);
    r_target = [-1;0;0]*0.95*Gravity.r_Ceres;

    Constraints.set_h(Constraints.h_land);
    Constraints.set_hdot(Constraints.hdot_land);
    Constraints.set_mu_phi(Gravity.mu_Ceres);
    Constraints.set_rho(Gravity.r_Ceres);
    Constraints.set_w_xmax(w_xmax);
    Constraints.set_w_umax(w_umax);
    Constraints.set_u_max(u_max);
    Constraints.set_obstacle(Gravity.obstacle_Ceres);
    Constraints.set_gravity(Gravity.gravity_Ceres);
    epsilon1 = -1/2*Constraints.PhiInv( delta + Constraints.Phi(0) - 2*(l_h*w_xmax + 1/2*gamma1)^2 );

    Disturbances.set_d_u(3);
    Disturbances.set_d_x(3);
    Disturbances.read_data();
    Disturbances.set_w_umax(w_umax);
    Disturbances.set_w_xmax(w_xmax);
    Disturbances.set_t_inc(10);

    Controllers.set_u_max(u_max);
    Controllers.set_w_umax(w_umax);
    Controllers.set_w_xmax(w_xmax);
    Controllers.set_r_target(r_target);
    Controllers.set_epsilon1(epsilon1);

    xdot_params = [(t,x)->Controllers.CalculateU(t,x,Constraint_func,case);
        Gravity.gravity_Ceres;
        Disturbances.Matched;
        Disturbances.Unmatched];
    ode_data = SavedValues(Float64, Vector{Float64});
    saving = SavingCallback(save_ode_data_landing, ode_data);
    stop_sim_func = (x,t,integrator)->Constraints.h(t,x) >= 0 || quit_early(x,t,integrator);
    stopping = DiscreteCallback(stop_sim_func,terminate!)
    problem = ODEProblem(UpdateX_Ceres!, x0, (0.0, 20000.0), xdot_params);
    sol = DifferentialEquations.solve(problem, reltol=1e-8, abstol=1e-7, isoutofdomain=DomainLimiter, dtmin=1e-10, callback=CallbackSet(saving,stopping));
elseif case==:Docking
    Controllers.set_f(Constraints.f_dock);
    Controllers.set_g(Constraints.g_dock);

    rT = [Gravity.R_Earth + Gravity.altitude_Earth; 0; 0];
    vT = [0; sqrt(Gravity.mu_Earth/rT[1]); 0];
    distance = 10e3;
    arc_length = 2*asin((distance/2) / rT[1]);
    rot = [cos(arc_length) sin(arc_length) 0; -sin(arc_length) cos(arc_length) 0; 0 0 1];
    rC = rot*rT;
    vC = rot*vT;

    rCT = rC - rT;
    omegaTN = [0;0;-Gravity.n];
    vCT = vC - vT + cross(omegaTN, rC-rT); # velocity as observed in Hill's rotating frame is zero

    x0 = [rCT[1:2]; vCT[1:2]];
    u_max = 0.082;
    w_umax = 0.002;
    w_xmax = 0.001;
    gamma1 = 0.07;
    gamma2 = 0.12;
    l_h = 1;
    delta = gamma2^2/2;
    if delta < 2*(l_h*w_xmax + gamma1/2)^2
        println("Warning: Landing cannot be achieved for these tolerances");
    end
    Delta = 0.03;
    v_max = 10;

    a_max_dock = u_max - 2*Gravity.n*v_max - w_umax;
    a_max_lr = u_max - 2*Gravity.n*v_max - 3*Gravity.n^2*distance;

    Constraint_func = (t,x)->Constraints.H_docking(t,x,delta);

    Constraints.set_Delta(Delta);
    Constraints.set_v_max(v_max);
    Constraints.set_w_xmax(w_xmax);
    Constraints.set_w_umax(w_umax);
    Constraints.set_u_max(u_max);
    Constraints.set_a_max_dock(a_max_dock);
    Constraints.set_a_max_lr(a_max_lr);
    Constraints.set_gravity(Gravity.HCW);
    alpha_inv_2_dock = 1/a_max_dock*( delta - 1/2*(2*l_h*w_xmax + gamma1)^2 );
    epsilon1 = alpha_inv_2_dock/2;
    epsilon_lr = 0.005; # This number must be sufficiently small that the QP is always feasible
    epsilon_v = 0.005*v_max; # arbitrarily set to 0.5% of the maximum v value

    Disturbances.set_d_u(2);
    Disturbances.set_d_x(2);
    Disturbances.read_data();
    Disturbances.set_w_umax(w_umax);
    Disturbances.set_w_xmax(w_xmax);
    Disturbances.set_t_inc(5);

    Controllers.set_u_max(u_max);
    Controllers.set_w_umax(w_umax);
    Controllers.set_w_xmax(w_xmax);
    Controllers.set_epsilon1(epsilon1);
    Controllers.set_epsilon_lr(epsilon_lr);
    Controllers.set_epsilon_v(epsilon_v);
    Controllers.set_enforce_h_left(DockingSwitch.enforce_h_left)

    xdot_params = [(t,x)->Controllers.CalculateU(t,x,Constraint_func,case);
        Gravity.HCW;
        Disturbances.Matched;
        Disturbances.Unmatched];
    ode_data = SavedValues(Float64, Vector{Float64});
    saving = SavingCallback(save_ode_data_docking, ode_data);
    Constraints.set_h(Constraints.h_dock);
    stop_sim_func = (x,t,integrator)->Constraints.h(t,x) >= 0 || quit_early(x,t,integrator);
    stopping = DiscreteCallback(stop_sim_func,terminate!)
    problem = ODEProblem(UpdateX_HCW!, x0, (0.0, 20000.0), xdot_params);
    sol = DifferentialEquations.solve(problem, reltol=1e-11, abstol=1e-10, isoutofdomain=DockingSwitch.DomainLimiter, dtmin=1e-10, callback=CallbackSet(saving,stopping));
end

function save_data(filename)
    file = open(filename, "w+");
    for i=1:length(ode_data.t)
        @printf(file, "%.8f ", ode_data.t[i])
        @printf(file, "%.8f %.8f %.8f %.8f %.8f %.8f ", ode_data.saveval[i][1], ode_data.saveval[i][2], ode_data.saveval[i][3], ode_data.saveval[i][4], ode_data.saveval[i][5], ode_data.saveval[i][6])
        @printf(file, "%.8f ", ode_data.saveval[i][7]) # H
        @printf(file, "%.8f %.8f %.8f ", ode_data.saveval[i][8], ode_data.saveval[i][9], ode_data.saveval[i][10]) # u
        @printf(file, "%.8f ", ode_data.saveval[i][11]) # sigma
        @printf(file, "%.8f ", ode_data.saveval[i][12]) # beta
        @printf(file, "%.8f\n", ode_data.saveval[i][13]) # residual
    end
    close(file);
end

function save_control_data(filename)
    file = open(filename, "w+");
    for i=1:length(ode_data.t)
        @printf(file, "%.8f ", ode_data.t[i])
        u = xdot_params[1](ode_data.t[i], ode_data.saveval[i][1:length(x0)])
        if case==:Landing
            @printf(file, "%.8f ", Constraints.last_H)
            @printf(file, "%.8f %.8f %.8f\n", u[1], u[2], u[3]) # real control inputs instead of "last control inputs"
        elseif case==:Docking
            @printf(file, "%.8f %.8f %.8f %.8f ", Constraints.last_H[1], Constraints.last_H[6], Constraints.last_H[7], Constraints.last_H[2])
            @printf(file, "%.8f %.8f\n", u[1], u[2]) 
        end
    end
    close(file);
end

save_data("Sim.txt");
save_control_data("Control.txt"); # this takes a while to run

function move_data(name)
    Base.Filesystem.cp("Sim.txt", "Results/Sim"*name*".txt", force=true);
    Base.Filesystem.cp("Control.txt", "Results/Control"*name*".txt", force=true);
end

# if case==:Landing
#     move_data("Landing");
# elseif case==:Docking
#     move_data("Docking");
# end

println("\007")
