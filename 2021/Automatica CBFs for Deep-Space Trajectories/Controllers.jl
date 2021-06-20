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
using_switching = true;

function set_u_max(x); global u_max = x; end
function set_w_umax(x); global w_umax = x; end
function set_w_xmax(x); global w_xmax = x; end
function set_using_switching(x::Bool); global using_switching = x; end

r_target = [0.0;0.0;0.0];
v_target = [0.0;0.0;0.0];
function set_r_target(x); global r_target = x; end
function set_v_target(x); global v_target = x; end

last_H = NaN;
last_u = [NaN; NaN; NaN];
last_sigma = NaN;
last_residual = NaN;

epsilon1 = NaN;
epsilon2 = NaN;
slope = NaN;
function set_epsilon1(x); global epsilon1 = x; end
function set_epsilon2(x); global epsilon2 = x; end
function set_slope(x); global slope = x; end

function CalculateU(t, x, Constraint_func, sim_case, switching_func)

function alpha(lambda, W=NaN)
    if using_switching
        # @show W/epsilon1
        global last_residual = W/epsilon1;
        return W*lambda/epsilon1;
    else
        return slope*lambda;
    end
end

# CBF Conditions
H, partial_H, grad_H = Constraint_func(t, x);
n = size(grad_H)[1];
A = zeros(n, 3);
b = zeros(n)
for i=1:n
    W = norm( grad_H[i,:]'*g(t,x) )*w_umax + norm( grad_H[i,:] )*w_xmax;
    b[i] = alpha(-H[i],W) - partial_H[i] - dot(grad_H[i,:], f(t,x)) - norm(grad_H[i,:]'*g(t,x))*w_umax - norm(grad_H[i,:])*w_xmax;
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

# Nominal controller
r = x[1:3];
v = x[4:6];

if sim_case==:Flyby
    kp = 0.00003*(u_max/0.5)^2/100*1e3;
    kd = 0.03*(u_max/0.5)/5*5*10
    v_target[1] = sqrt(2*Gravity.mu_Ceres/norm(r) + 100^2);
    global r_target = dot(r,v_target)/dot(v_target,v_target)*v_target;
elseif sim_case==:Landing
    kp = 0.0000075;
    kd = 0.03;
elseif sim_case==:Eros
    kp = 0.00003;
    kd = 0.03;
end
u_nom = kp*(r_target - r) + kd*(v_target - v);
F = -2*u_nom;
J = [1.0 0.0 0.0;
     0.0 1.0 0.0;
     0.0 0.0 1.0];
u = Convex.Variable(3)
lower = -[u_max; u_max; u_max];
upper =  [u_max; u_max; u_max];

if sim_case == :Eros
    sliding = n >= 1 || !using_switching
    if n >= 1
        global last_H = maximum(H);
    else
        global last_H = -Base.Inf;
    end
else
    sliding = switching_func(t,H)[1] || !using_switching;
    global last_H = H;
end

if sim_case == :Flyby
    scale = 10/u_max;
    A *= 1e6;
    b *= 1e6*scale;
    lower *= scale;
    upper *= scale;
    F *= scale;
end

if sliding
    problem = minimize(quadform(u, J) + dot(F,u), A*u <= b, u[1:3] <= upper, u[1:3] >= lower);
    global last_sigma = 1.0;
else
    problem = minimize(quadform(u, J) + dot(F,u), u[1:3] <= upper, u[1:3] >= lower);
    global last_sigma = 0.0;
end
Convex.solve!(problem, () -> ECOS.Optimizer(), verbose=false, silent_solver=true);

@show t, last_H

if n > 0 && sliding && maximum(A*u.value - b) > 1e-6
    # @show maximum(A*u.value - b)
    # Check if the quadratic program solver is working correctly or not
end

# If the QP solver fails and produces a very large control input, it is sometimes necessary
# to limit that control input so the ODE propagator works properly.
# function limit(x)
#     if x > 11
#         return 11;
#     elseif x < -11
#         return -11;
#     else
#         return x;
#     end
# end
# if sim_case == :Flyby && Main.cbf == :Integrated
#     u.value = limit.(u.value);
# end

if sim_case == :Flyby
    global last_u = u.value/scale;
    return u.value[1:3]/scale
else
    global last_u = u.value
    return u.value[1:3]
end

end # end function

end # end module

module Switching

time_start_sliding = [];
time_end_sliding = [];
sliding = [];

max_count = 0;

function switching_func(t,H,epsilon1,epsilon2)
    if t==0.0
        global time_start_sliding = -2*ones(length(H));
        global time_end_sliding = -1*ones(length(H));
        global sliding = Bool.(ones(length(H)));
    end
    count = 0;
    for i=1:length(H)
        # The control law for sigma is more complicated than in the paper because the 
        # variable step integrator could step backwards in time.
        if H[i] >= -epsilon1
            global sliding[i] = true;
            count += 1;
            global time_start_sliding[i] = t;
            global time_end_sliding[i] = Base.Inf;
        elseif H[i] < -epsilon2
            global sliding[i] = false;
            global time_end_sliding[i] = t;
        else
            if time_start_sliding[i] < time_end_sliding[i] # if both of these variables have already been set and we've back-stepped somewhere in between them
                if t >= time_start_sliding[i] && t < time_end_sliding[i]
                    global sliding[i] = true;
                    count += 1;
                else
                    global sliding[i] = false;
                end
            else # if we have not yet identified when sliding should stop
                if t >= time_start_sliding[i] # keep sliding if we already started sliding
                    global sliding[i] = true;
                    count += 1;
                else # stop sliding if we back-stepped before sliding started
                    global sliding[i] = false;
                    # this condition is why time_end_sliding > time_start_sliding when initialized
                end
            end
        end
    end
    global max_count = max(count, max_count)
    # println(max_count);
    return sliding;
end

end

nothing