# A file for all the constraint functions that could be needed

module Constraints

using LinearAlgebra
using ForwardDiff
using DifferentialEquations
using SpecialMatrices
include("Utilities.jl")
include("Gravity.jl")

mu = NaN;
u_max = NaN;
rho = NaN;
w_umax = NaN;
w_xmax = NaN;
a_max = NaN;

last_H = NaN;
last_t = NaN;

Z3 = [0.0 0.0 0.0;
      0.0 0.0 0.0;
      0.0 0.0 0.0];
I3 = [1.0 0.0 0.0;
      0.0 1.0 0.0;
      0.0 0.0 1.0];

obstacle = nothing;
gravity = nothing;
partial_gravity = nothing;
grad_gravity = nothing;

function Phi(h)
    return mu/(rho - h) + (w_umax - u_max)*h
end

function PhiInv(y)
    c1 = u_max - w_umax
    Phi_crit = 2*sqrt(mu*c1) - c1*rho
    # if y >= Phi_crit
    try
        return (c1*rho - y - sqrt((y - c1*rho)^2 - 4*c1*(mu - rho*y)))/(2*c1)
    # else
    catch
        return rho - sqrt(mu/c1)+10 # to deal with potentially undefined values of PhiInv
    end
end

function phi(h)
    return mu/(rho - h)^2 + w_umax - u_max
end

###################################################################################################
###################################################################################################

function h(t,x,args=[])
    r = x[1:3];
    v = x[4:6];
    rc = obstacle(t,args)[1];
    return rho - norm(r-rc);
end

function hdot(t,x,args=[])
    r = x[1:3];
    v = x[4:6];
    rc, vc = obstacle(t,args)[1:2];
    return -dot(r-rc, v-vc)/norm(r-rc);
end

function partial_h(t,x,args=[])
    r = x[1:3];
    v = x[4:6];
    rc, vc = obstacle(t,args)[1:2];
    return -dot(r-rc, -vc)/norm(r-rc);
end

function grad_h(t,x,args=[])
    r = x[1:3];
    v = x[4:6];
    rc = obstacle(t,args)[1];
    return [-(r-rc)'/norm(r-rc) [0.0 0.0 0.0]];
end

function hddot(t,x,u,args=[])
    r = x[1:3];
    v = x[4:6];
    rc, vc, uc = obstacle(t,args)[1:3];
    g = gravity(t, r);
    return -norm(cross(r-rc,v-vc))^2/norm(r)^3 - dot(r-rc, u-uc+g)/norm(r)
end

function hdddot(t,x,u,udot,args=[])
    r = x[1:3];
    v = x[4:6];
    rc, vc, uc, udotc = obstacle(t,args);
    g = gravity(t, r);
    gdot = partial_gravity(t, r) + [grad_gravity(t, r) Z3]*f(t,x); # By construction, grad_gravity_central*g=0, so we omit that computation
    return 2*dot(v-vc, skew(r-vc)^2*(u-uc+g))/norm(r-rc)^3 + 3*norm(cross(r-rc,v-vc))^2*dot(r-rc,v-vc)/norm(r-rc)^5 - 
        dot(v-vc,u-uc+g)/norm(r-rc) + dot(r-rc,v-vc)*dot(r-rc,u-uc+g)/norm(r-rc)^3 - dot(r-rc,udot-udotc+gdot)/norm(r-rc);
end

###################################################################################################
###################################################################################################

function f(t,x)
    r = x[1:3];
    v = x[4:6];
    return [v; gravity(t, r)];
end

function g(t,x)
    return [0.0 0.0 0.0;
            0.0 0.0 0.0;
            0.0 0.0 0.0;
            1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0];
end

###################################################################################################
###################################################################################################

function odefunc!(xdot, x, controller, t)
    r = x[1:3];
    v = x[4:6];
    u = controller(t, x);
    vdot = u + gravity(t,r);
    xdot[1:3] = v;
    xdot[4:6] = vdot;
end

function grad_h_odefunc!(xdot, x, controller, t)
    y = x[1:6];
    Theta = reshape(x[7:42], (6, 6));
    ydot = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
    odefunc!(ydot, y, controller, t);
    df_dy = [Z3 I3; grad_gravity(t,x[1:3]) Z3];
    df_dt = [0.0; 0.0; 0.0; partial_gravity(t,x[1:3])];
    # g(t,x) is a constant, so dg_dy = dg_dt = 0
    du_dy = ForwardDiff.jacobian(z -> controller(t,z), y);
    du_dt = ForwardDiff.derivative(s -> controller(s,x), t);
    Thetadot = (df_dy + g(t,x)*du_dy)*Theta;
    thetadot = df_dt + g(t,x)*du_dt;
    xdot[1:6] = ydot;
    xdot[7:42] = reshape(Thetadot, (36, 1));
    xdot[43:48] = thetadot;
end
# A test run had the exact jacobian of gravity_central taking 1 us to run, and forward diff taking 7 us to run. 
# Thus, it is slightly better to keep the derivatives explicit, but for simplicity, I'm only doing that with certain functions

function ode_stop_event3(x, t, integrator, controller)
    u = controller(t,x);
    udot = ForwardDiff.derivative(s -> controller(s,x), t) + ForwardDiff.jacobian(y -> controller(t,y), x)*(f(t,x)+g(t,x)*u);
    c3 = hdddot(t,x,u,udot);
    if c3 > 0
        value = c3;
    else
        c2 = hddot(t,x,u);
        if c2 > 0
            value = c2;
        else
            value = hdot(t,x) + 100; # ensure there is enough to curve fit
        end
    end
    return value < 0;
end

function ode_stop_event2(x, t, integrator, controller)
    u = controller(t,x);
    c2 = hddot(t,x,u);
    if c2 > 0
        value = c2;
    else
        value = hdot(t,x) + 100; # ensure there is enough to curve fit
    end
    return value < 0;
end

###################################################################################################
###################################################################################################

function H_two_norm(t,x)
    function func(t,x)
        return h(t,x) + absSq(hdot(t,x) + w_xmax)/(2*a_max);
    end
    H = func(t,x);
    partial_H = ForwardDiff.derivative(y -> func(y,x), t);
    grad_H = ForwardDiff.gradient(y -> func(t,y), x)';
    global last_H = H;
    return H, partial_H, grad_H;
end

function H_two_norm_mesh(t,x,switching_func)
    function func(t,x,i)
        return h(t,x,i) + absSq(hdot(t,x,i) + w_xmax)/(2*a_max);
    end
    H = zeros(Gravity.Eros_num_points);
    for i=1:Gravity.Eros_num_points
        H[i] = func(t,x,i);
    end
    switch_cond = switching_func(H)
    indices = findall(switch_cond);
    H_select = H[indices];
    partial_H = zeros(length(indices));
    grad_H = zeros(length(indices), 6);
    for i=1:length(indices)
        partial_H[i] = ForwardDiff.derivative(y -> func(y,x,indices[i]), t);
        grad_H[i,:] = ForwardDiff.gradient(y -> func(t,y,indices[i]), x)';
    end
    global last_H = maximum(H);
    return H_select, partial_H, grad_H 
end

function H_variable(t,x;delta=0)
    function func(t,x)
        return PhiInv( Phi(h(t,x)) - absSq(hdot(t,x) + w_xmax)/2 + delta );
    end
    H = func(t,x);
    partial_H = ForwardDiff.derivative(y -> func(y,x), t);
    grad_H = ForwardDiff.gradient(y -> func(t,y), x)';
    global last_H = H;
    return H, partial_H, grad_H;
end


"""
4 milliseconds to run this function
"""
function H_predictive(t,x,p)
    controller = p[1];
    stop_func = p[2];
    duration = p[3];
    problem = ODEProblem(odefunc!, x, (0.0, duration), controller);
    end_condition = DiscreteCallback(stop_func,terminate!);
    # if duration > 5e5
        sol = solve(problem, reltol=1e-12, abstol=1e-8, callback=end_condition);
    # else
    #     if last_t >= 500
    #         sol = solve(problem, reltol=1e-16, abstol=1e-10, callback=end_condition);
    #     else
    #         sol = solve(problem, reltol=1e-17, abstol=1e-11, callback=end_condition);
    #     end
    # end
    hvals = h.(sol.t, sol.u);

    if hvals[2] < hvals[1]
        d = diff(hvals) .> 0;
        if sum(d) >= 1
            j = findmax(d)[2];
        else
            j = 1;
        end
    else
        j = 1;
    end

    H, i = findmax(hvals[j:end]);
    i = i+j-1;
    t_end = sol.t[i];
    if t_end > 0.0
        if i==1; i = 2; end
        if i==length(sol.t); i=i-1; end
        t1 = sol.t[i-1];
        t2 = sol.t[i];
        t3 = sol.t[i+1];
        y1 = hvals[i-1];
        y2 = hvals[i];
        y3 = hvals[i+1];
        M = Vandermonde([t1; t2; t3]);
        yarr = [y1;y2;y3];
        cba = \(M, yarr); # quadratic polynomial coefficients
        if cba[3] < 0
            t_end = max(0, -cba[2]/(2*cba[3]));
            H = cba[3]*t_end^2 + cba[2]*t_end + cba[1];
        end
    end
    global last_t = t_end
    if t_end >= duration
        println("t_end = duration. Solution likely innaccurate.");
    end

    # Next step is to find the derivatives
    I6 = [I3 Z3; Z3 I3];
    x0 = [x; reshape(I6, (36, 1)); 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]; # 48-dimensional: [y; Theta; theta]
    if H <= 0
        if t_end==0
            y = x0;
        elseif t_end > 0
            problem = ODEProblem(grad_h_odefunc!, x0, (0.0, t_end), controller);
            if duration > 5e5
                sol = solve(problem, reltol=1e-6);
            else
                sol = solve(problem, reltol=1e-8); # the really tight tolerance is not needed here
            end
            y = last(sol.u);
        else
            error("Invalid t. Check H");
        end

        # The following formulas are true in either the zero or nonzero case
        x_end = y[1:6];
        grad_H = grad_h(t_end, x_end) * reshape(y[7:42], (6, 6));
        partial_H = partial_h(t_end, x_end) + dot(grad_h(t_end, x_end), y[43:48]);
    else
        partial_H = 0.0;
        grad_H = [1.0 1.0 1.0 1.0 1.0 1.0];
        # This means we're going to reject this time step anyways, so we might as well save integration time
    end

    global last_H = H;
    return H, partial_H, grad_H
end

###################################################################################################
###################################################################################################

function set_mu_phi(x); global mu = x; end
function set_rho(x); global rho = x; end
function set_u_max(x); global u_max = x; end
function set_a_max(x); global a_max = x; end
function set_w_umax(x); global w_umax = x;end
function set_w_xmax(x); global w_xmax = x; end
function set_gravity(x); global gravity = x; end
function set_partial_gravity(x); global partial_gravity = x; end
function set_grad_gravity(x); global grad_gravity = x; end
function set_obstacle(x); global obstacle = x; end

function controller2(t, x)
    rc = obstacle(t, [])[1];
    r = x[1:3]
    return (u_max - w_umax)*(r - rc)/norm(r - rc);
end

function controller2S(t, x)
    rc = obstacle(t, [])[1];
    r = x[1:3]
    return (u_max - w_umax)*(r - rc)/norm(r - rc, Inf);
end

function controller3(t, x)
    rc, vc = obstacle(t, [])[1:2];
    r = x[1:3];
    v = x[4:6];
    v_orth = v-vc - dot(r-rc, v-vc)/dot(r-rc, r-rc)*(r-rc);
    return (u_max - w_umax)*v_orth/(norm(v_orth)+1e-8);
end

function controller3S(t, x)
    rc, vc = obstacle(t,[])[1:2];
    r = x[1:3];
    v = x[4:6];
    v_orth = v-vc - dot(r-rc, v-vc)/dot(r-rc, r-rc)*(r-rc);
    return (u_max - w_umax)*v_orth/(norm(v_orth, Inf)+1e-8);
end

end

nothing