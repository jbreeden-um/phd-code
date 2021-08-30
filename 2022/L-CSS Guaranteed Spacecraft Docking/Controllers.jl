module Controllers

include("Gravity.jl");

f = nothing;
g = nothing;
function set_f(x); global f = x; end
function set_g(x); global g = x; end

using Convex
using ECOS

u_max = NaN;
w_umax = NaN;
w_xmax = NaN;

function set_u_max(x); global u_max = x; end
function set_w_umax(x); global w_umax = x; end
function set_w_xmax(x); global w_xmax = x; end

r_target = [0.0;0.0;0.0];
v_target = [0.0;0.0;0.0];
function set_r_target(x); global r_target = x; end
function set_v_target(x); global v_target = x; end

last_H = NaN;
last_H_left = NaN;
last_H_right = NaN;
last_u = [NaN; NaN; NaN];
last_sigma = NaN;
last_residual = NaN;

enforce_h_left = nothing;
function set_enforce_h_left(x); global enforce_h_left = x; end

epsilon1 = NaN;
epsilon_lr = NaN;
epsilon_v = NaN;
function set_epsilon1(x); global epsilon1 = x; end
function set_epsilon_lr(x); global epsilon_lr = x; end
function set_epsilon_v(x); global epsilon_v = x; end

function CalculateU(t, x, Constraint_func, sim_case)

function alpha(lambda, W=NaN)
    return W*lambda/epsilon1;
end

# CBF Conditions
H, partial_H, grad_H = Constraint_func(t, x);
global last_H = H[1];
n = 0;
m = 0;
if sim_case == :Landing
    n = 1;
    m = 3;
elseif sim_case == :Docking
    global last_H_right = H[6];
    global last_H_left = H[7];
    H[2:5] *= epsilon1 / epsilon_v; # This is for the alpha function having a different slope
    H[6] *= epsilon1 / epsilon_lr; 
    H[7] *= epsilon1 / epsilon_lr;
    if enforce_h_left(t)
        n = 7;
    else
        n = 6;
    end
    m = 2;
end
A = zeros(n, m);
b = zeros(n)
for i=1:n
    # The last term of the following line only considers the 1:3 indices because that is how we defined our specific model.
    # The separation of w_x and w_u into different rows makes sense because norm(grad_H[4:6]) >> norm(grad_H[1:3]).
    # If you are copying this code for other purposes, to follow the general formulation in the paper, replace grad_H[i,1:3] with grad_H[i,:].
    global W = norm( grad_H[i,:]'*g(t,x) )*w_umax + norm( grad_H[i,1:m] )*w_xmax;

    b[i] = alpha(-H[i],W) - partial_H[i] - dot(grad_H[i,:], f(t,x)) - W;
    A[i,:] = grad_H[i,:]'*g(t,x);

    if norm(A[i,:]) == 0 && b[i] < 0
        println("LgH is zero, but the constraint is still active.");
        @show t
        @show x
        @show A
        @show b
    end

    scale = eps()/maximum(abs.(A[i,:]))*1e12;
    if !isinf(scale)
        A[i,:] *= scale;
        b[i] *= scale;
    end
end

if sim_case == :Landing
    u_nom = A[1:3]/norm(A)^2;

    c = Convex.Variable(1)
    upper = u_max / maximum(abs.(u_nom))
    lower = -upper

    problem = maximize(c, A*u_nom*c <= b, c[1] <= upper, c[1] >= lower);
    
    Convex.solve!(problem, () -> ECOS.Optimizer(), verbose=false, silent_solver=true);
    u = u_nom*c.value[1];
    
    if maximum(abs.(u)) - u_max > 1e-5 && last_H <= 0
        println("Something is wrong")
        @show b[1], c.value[1], u
    end

    # The following line works great, but uses large control inputs for the first quarter second
    # u = u_nom;

    @show t, last_H
    if maximum(A*u-b) > 1e-10 && last_H <= 0
        @show b, A*u-b
    end

    function limit(x, m)
        if x > m
            return m;
        elseif x < -m
            return -m;
        else
            return x;
        end
    end

    if last_H > 0
        # The time step will be discarded anyways, but we want to make sure u is reasonably conditioned because of the Runge-Kutta integrator.
        # Note that we only apply this check when last_H > 0. If u is too large for safe H values, then we want to know about it, so we don't fix it here.
        u = limit.(u, u_max)
    end

    global last_u = u;
    return u;
elseif sim_case == :Docking
    u_nom1 = b[1]*A[1,:]/norm(A[1,:])^2;
    kp = 0.1;
    u_nom2 = [-kp*x[1]; 0];
    # u_nom1 and u_nom2 are always orthogonal so we don't have to worry as much about Lipschitz continuity of the QP as with a pure linear control law.

    # @show u_nom1, u_nom2
    F = -2*(u_nom1 + u_nom2);
    J = [1.0 0.0;
        0.0 1.0];
    u = Convex.Variable(2)
    lower = -[u_max; u_max];
    upper =  [u_max; u_max];
    problem = minimize(quadform(u, J) + dot(F,u), A*u <= b, u <= upper, u >= lower);
    Convex.solve!(problem, () -> ECOS.Optimizer(), verbose=false, silent_solver=true);

    if Int(problem.status) == 2 && last_H <= 0
        println("Warning: Optimization Error: ")
        for i=1:n
            print("|A| = "*string(norm(grad_H[i,:]'*g(t,x))))
            W = norm( grad_H[i,:]'*g(t,x) )*w_umax + norm( grad_H[i,1:m] )*w_xmax;
            print(" . alpha = "*string(alpha(-H[i],W)))
            print(" . Lf = "*string(- dot(grad_H[i,:], f(t,x))))
            println(" . W = "*string(-W))
        end
    end

    # @show A*u.value - b

    if maximum(A*u.value - b) > 5e-6
        @show maximum(A*u.value - b)
        # Check if the quadratic program solver is working correctly or not
    end
    @show t, last_H
    # @show t, last_H, u.value, n


    global last_u = u.value;
    return u.value;
end

end # end function

end # end module

"""
This module is used to determine the appropriate time to start enforcing both the left side
and right side constraints.
    enforce_h_left() returns whether the controller should consider both constraints
    DomainLimiter() returns whether the current step is valid by considering whether both constraints need to be enforced or not.
"""
module DockingSwitch
time_enter = Base.Inf;

function enforce_h_left(t)
    if Main.Controllers.last_H_left <= 0 && Main.Controllers.last_H_right <=0 && t < time_enter
        global time_enter = t;
    end
    if t < time_enter
        return false
    else
        return true
    end
end

function DomainLimiter(x, p, t)
    if Main.Controllers.last_H_left <= 0 && Main.Controllers.last_H_right <=0 && t < time_enter
        global time_enter = t;
    end
    con = false;
    # The right side constraint (too large of an orbit radius) is always enforced, while the left side constraint (too small of an orbit radius) is only enforced once initial alignment with the docking axis is complete.
    if t < time_enter
        con = Main.Controllers.last_H_right <=0;
    else
        con = Main.Controllers.last_H_left <= 0 && Main.Controllers.last_H_right <=0;
    end
    return Main.Controllers.last_H > 0 || isnan(Main.Controllers.last_u[1]) || !con;
end
end

nothing