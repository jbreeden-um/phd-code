
using DifferentialEquations
using Plots
using Printf
using LinearAlgebra

try
    Controllers
catch
    include("Dynamics.jl"); include("Constraints.jl"); include("Controllers.jl"); include("Disturbances.jl");
end

try 
    if !(AutonomousExecution)
        error("Using manual settings.")
    end
catch
    # CHANGE THESE SETTINGS TO SWITCH BETWEEN SIMULATION CASES
    # global case = :Eros
    # global case = :Landing
    global case = :Flyby
        # global cbf = :Constant
        # global cbf = :Variable
        # global cbf = :Integrated
        # global cbf = :Orbital
        # global cbf = :NewAlpha
        # global cbf = :NoSwitch
        global cbf = :None
end

# The reader is welcome to test this code with different parameters, settings, tolerances, etc.
# The code below describes what was done to generate the plots in the paper.
# There are also a few commented out sections of code, which should be treated as tools
# that may improve execution when settings are changed.

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

function save_ode_data(x,t,integrator)
    return [x; Controllers.last_H; Controllers.last_u[1:3]; Controllers.last_sigma; Constraints.last_t; Controllers.last_residual];
end

Controllers.set_f(Constraints.f);
Controllers.set_g(Constraints.g);

if case==:Flyby
    r0 = [-6e7; -1e6; 0.0];
    v0 = [20.0; -2.0; 0.0];
    x0 = [r0; v0];

    u_max = 1e-4; # approximate acceleration of an ion drive
        # 91 mN at 982 kg https://solarsystem.nasa.gov/missions/dawn/technology/ion-propulsion/
        # 1.15 mN at 12 kg http://www.busek.com/technologies__ion.htm
    w_umax = 5e-6;
    w_xmax = 2e-6;

    epsilon1 = 5e4;
    epsilon2 = 3*epsilon1;
    v_target = [NaN;0.0;0.0];

    Constraints.set_mu_phi(Gravity.mu_Ceres);
    if cbf==:Constant
        Constraint_func = Constraints.H_two_norm;
        r_min = sqrt(Gravity.mu_Ceres/(0.5*(u_max-w_umax))); # r_min was set arbitrarily as the place at which a_max = 0.5 u_max
            # Note that a_max does not exist at Gravity.r_Ceres
        Constraints.set_rho(r_min);
        Constraints.set_a_max(0.5*u_max - 0.5*w_umax - w_xmax);
        rTol = 1e-5;
        aTol = 1e-4;
        Controllers.set_using_switching(true);
    elseif cbf==:NewAlpha
        Constraint_func = Constraints.H_two_norm;
        r_min = sqrt(Gravity.mu_Ceres/(0.5*(u_max-w_umax)));
        Constraints.set_rho(r_min);
        Constraints.set_a_max(0.5*u_max - 0.5*w_umax - w_xmax);
        rTol = 1e-5;
        aTol = 1e-4;
        Controllers.set_using_switching(false);
        Controllers.set_slope(x->x/epsilon1);
    elseif cbf==:NoSwitch
        Constraint_func = Constraints.H_two_norm;
        r_min = sqrt(Gravity.mu_Ceres/(0.5*(u_max-w_umax)));
        Constraints.set_rho(r_min);
        Constraints.set_a_max(0.5*u_max - 0.5*w_umax - w_xmax);
        rTol = 1e-5;
        aTol = 1e-4;
        Controllers.set_using_switching(false);
        # Controllers.set_slope(x->1.90e-4); # average during the nominal simulation
        Controllers.set_slope(x->3.42e-5); # average during the nominal simulation when sigma was active
        # Controllers.set_slope(x->1e-4); # largest slope I've tested that does not fail due to numerical tolerances, though the control input is bad
    elseif cbf==:Variable
        Constraint_func = Constraints.H_variable;
        r_min = 1.25*sqrt(Gravity.mu_Ceres/(u_max - w_umax)) # to ensure Phi does not become positive
        Constraints.set_rho(r_min);
        rTol = 1e-6;
        aTol = 1e-5;
        Controllers.set_using_switching(true);
    elseif cbf==:Integrated
        control = Constraints.controller2S;
        stop = (a,b,c)->Constraints.ode_stop_event2(a,b,c,control);
        param = (control, stop, 16e6)
        Constraint_func = (t,x)->Constraints.H_predictive(t,x,param);
        w_xmax = 0;
        Constraints.set_rho(Gravity.r_Ceres);
        r_min = sqrt(Gravity.mu_Ceres/(u_max - w_umax))
        Constraints.set_rho(r_min);
        rTol = 1e-4;
        aTol = 1e-3;
        Controllers.set_using_switching(true);
    elseif cbf==:Orbital
        control = Constraints.controller3S;
        stop = (a,b,c)->Constraints.ode_stop_event3(a,b,c,control);
        param = (control, stop, 8e6);
        Constraint_func = (t,x)->Constraints.H_predictive(t,x,param);
        w_xmax = 0;
        Constraints.set_rho(Gravity.r_Ceres);
        rTol = 1e-6;
        aTol = 1e-5;
        Controllers.set_using_switching(true);
    elseif cbf==:None
        Constraint_func = (t,x)->(-1, [], []);
        rTol = 1e-5;
        aTol = 1e-4;
        Controllers.set_using_switching(true);
    else
        error("Unknown CBF method")
    end
    
    Constraints.set_w_xmax(w_xmax);
    Constraints.set_w_umax(w_umax);
    Constraints.set_u_max(u_max);
    Constraints.set_obstacle(Gravity.obstacle_Ceres);
    Constraints.set_gravity(Gravity.gravity_Ceres);
    Constraints.set_partial_gravity(Gravity.partial_gravity_Ceres);
    Constraints.set_grad_gravity(Gravity.grad_gravity_Ceres);

    Disturbances.read_data();
    Disturbances.set_w_umax(w_umax);
    Disturbances.set_w_xmax(w_xmax);
    Disturbances.set_t_inc(10000);

    Controllers.set_u_max(u_max);
    Controllers.set_w_umax(w_umax);
    Controllers.set_w_xmax(w_xmax);
    Controllers.set_v_target(v_target);
    Controllers.set_epsilon1(epsilon1);
    Controllers.set_epsilon2(epsilon2);

    if cbf==:None
        switch_func = (s,y)->false
    else
        switch_func = (s,y)->Switching.switching_func(s,y,epsilon1,epsilon2)
    end
    xdot_params = [(t,x)->Controllers.CalculateU(t,x,Constraint_func,case,switch_func);
        Gravity.gravity_Ceres;
        Disturbances.Matched;
        Disturbances.Unmatched];
    ode_data = SavedValues(Float64, Vector{Float64});
    saving = SavingCallback(save_ode_data, ode_data);
    stopping = DiscreteCallback(quit_early,terminate!)
    if cbf == :Orbital || cbf==:Integrated || cbf==:NoSwitch
        problem = ODEProblem(UpdateX!, x0, (0.0, 6000000.0), xdot_params);
        sol = DifferentialEquations.solve(problem, reltol=rTol, abstol=aTol, isoutofdomain=DomainLimiter, dtmin=1e-8, callback=CallbackSet(saving,stopping));
    else
        problem = ODEProblem(UpdateX!, x0, (0.0, 6000000.0), xdot_params);
        sol = DifferentialEquations.solve(problem,ImplicitEuler(autodiff=false), isoutofdomain=DomainLimiter, reltol=rTol, abstol=aTol, dtmin=1e-10, callback=CallbackSet(saving,stopping));
    end
    r_target = [-x0[1]; 0; 0]; # for the plotting routine
elseif case==:Landing
    # "Landing" is no longer covered in this paper due to length constraints, but the code is still included.
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

    Constraints.set_mu_phi(Gravity.mu_Ceres);
    Constraints.set_rho(Gravity.r_Ceres);
    Constraints.set_w_xmax(w_xmax);
    Constraints.set_w_umax(w_umax);
    Constraints.set_u_max(u_max);
    Constraints.set_obstacle(Gravity.obstacle_Ceres);
    Constraints.set_gravity(Gravity.gravity_Ceres);
    epsilon1 = -1/2*Constraints.PhiInv( delta + Constraints.Phi(0) - 2*(l_h*w_xmax + 1/2*gamma1)^2 );
    epsilon2 = 3*epsilon1;

    Disturbances.read_data();
    Disturbances.set_w_umax(w_umax);
    Disturbances.set_w_xmax(w_xmax);
    Disturbances.set_t_inc(10);

    Controllers.set_u_max(u_max);
    Controllers.set_w_umax(w_umax);
    Controllers.set_w_xmax(w_xmax);
    Controllers.set_r_target(r_target);
    Controllers.set_epsilon1(epsilon1);
    Controllers.set_epsilon2(epsilon2);
    Controllers.set_using_switching(true);

    xdot_params = [(t,x)->Controllers.CalculateU(t,x,Constraint_func,case,(s,y)->Switching.switching_func(s,y,epsilon1,epsilon2));
        Gravity.gravity_Ceres;
        Disturbances.Matched;
        Disturbances.Unmatched];
    ode_data = SavedValues(Float64, Vector{Float64});
    saving = SavingCallback(save_ode_data, ode_data);
    stop_sim_func = (x,t,integrator)->Constraints.h(t,x) >= 0 || quit_early(x,t,integrator);
    stopping = DiscreteCallback(stop_sim_func,terminate!)
    problem = ODEProblem(UpdateX!, x0, (0.0, 20000.0), xdot_params);
    sol = DifferentialEquations.solve(problem, reltol=1e-8, isoutofdomain=DomainLimiter, dtmin=1e-10, callback=CallbackSet(saving,stopping));
elseif case==:Eros
    func = x->norm( Gravity.gravity_Eros(0,x) - cross(Gravity.Eros_omega_vec, cross(Gravity.Eros_omega_vec, x)) )
    a_asteroid = maximum(func.(Gravity.r_Eros));
    a_asteroid = ceil(a_asteroid*1e4)/1e4; # rounding for safety and numerics
    u_max = 0.1;
    w_umax = 0.005;
    w_xmax = 0.001;
    r_target = [20e3;0;0];
    x0 = [-20e3;-4000;0;1;1;0];
    epsilon1 = 100;
    epsilon2 = 300;
    Constraint_func = (t,x)->Constraints.H_two_norm_mesh(t,x,(h)->Switching.switching_func(t,h,epsilon1,epsilon2));

    Constraints.set_rho(500); # based off the distance between the mesh points
        # if this is too small relative to the mesh size, then the spacecraft could go between mesh points
        # maximum distance between the 3897 point mesh is 849.45 m. rho was set as half of this rounded up.
    # Constraints.set_rho(1000);
        # this choice of rho is intended to produce slightly smoother results; it sort of succeeds
    Constraints.set_w_xmax(w_xmax);
    Constraints.set_w_umax(w_umax);
    Constraints.set_a_max(u_max - a_asteroid - w_umax - w_xmax);
    Constraints.set_obstacle(Gravity.obstacle_Eros);
    Constraints.set_gravity(Gravity.gravity_Eros);

    Disturbances.read_data();
    Disturbances.set_w_umax(w_umax);
    Disturbances.set_w_xmax(w_xmax);
    Disturbances.set_t_inc(10);
    
    Controllers.set_u_max(u_max);
    Controllers.set_w_umax(w_umax);
    Controllers.set_w_xmax(w_xmax);
    Controllers.set_r_target(r_target);
    Controllers.set_epsilon1(epsilon1);
    Controllers.set_epsilon2(epsilon2);
    Controllers.set_using_switching(true);

    xdot_params = [(t,x)->Controllers.CalculateU(t,x,Constraint_func,case,nothing); Gravity.gravity_Eros; Disturbances.Matched; Disturbances.Unmatched];
    ode_data = SavedValues(Float64, Vector{Float64});
    saving = SavingCallback(save_ode_data, ode_data);
    stopping = DiscreteCallback(quit_early,terminate!)
    problem = ODEProblem(UpdateX!, x0, (0.0, 6000.0), xdot_params);
    sol = DifferentialEquations.solve(problem, reltol=1e-5, isoutofdomain=DomainLimiter, dtmin=1e-10, callback=CallbackSet(saving,stopping));
end

if cbf != :None
    Switching.end_sim(); # this closes the Switching.txt file so it can be overwritten by save_control_data()
end

function extract_dimension(data, d)
    return data[d];
end

x = extract_dimension.(sol.u, 1);
y = extract_dimension.(sol.u, 2);
z = extract_dimension.(sol.u, 3);
xplt = plot3d(x, y, z, reuse = false, xlabel="X", ylabel="Y", zlabel="Z", camera=(0,90), legend=:topleft);
display(xplt)

function closest(x,y,z)
    n = sqrt(x^2 + y^2 + z^2);
    out = [x y z]*Gravity.r_Ceres/n;
end
circ = closest.(x,y,z)
c1 = extract_dimension.(circ, 1);
c2 = extract_dimension.(circ, 2);
c3 = extract_dimension.(circ, 3);
plot!(xplt, c1, c2, c3)
plot!(xplt, [r_target[1]], [r_target[2]], [r_target[3]], markershape=:circle, markersize=6, linealpha=0);
gui()

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
        u = xdot_params[1](ode_data.t[i], ode_data.saveval[i][1:6])
        @printf(file, "%.8f ", Constraints.last_H)
        @printf(file, "%.8f %.8f %.8f\n", u[1], u[2], u[3]) # real control inputs instead of "last control inputs"
    end
    close(file);
    if cbf != :None
        Switching.end_sim();
    end
end

save_data("Sim.txt");
save_control_data("Control.txt"); # this takes a while to run

println("\007")
