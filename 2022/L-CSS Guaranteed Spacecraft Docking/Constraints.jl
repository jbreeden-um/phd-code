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
Delta = NaN;
w_umax = NaN;
w_xmax = NaN;
a_max_dock = NaN;
a_max_lr = NaN;

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
v_max = nothing;

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

h = nothing;
hdot = nothing;

function h_land(t,x,args=[])
    r = x[1:3];
    v = x[4:6];
    rc = obstacle(t,args)[1];
    return rho - norm(r-rc);
end

function hdot_land(t,x,args=[])
    r = x[1:3];
    v = x[4:6];
    rc, vc = obstacle(t,args)[1:2];
    return -dot(r-rc, v-vc)/norm(r-rc);
end

# Start behind the target in its orbit, so the tangential coordinate is negative
function h_dock(t,x)
    return x[2];
end

function hdot_dock(t,x)
    return x[4];
end

# Start behind the target in its orbit, so the radial coordinate is negative. This constraint will initially be violated.
# Left as viewed from the chaser looking towards the target where orbit normal is the up direction (i.e. left = negative e_r direction in e_r/e_n plane)
function h_dock_left(t,x)
    return -x[1] - Delta;
end

function hdot_dock_left(t,x)
    return -x[3];
end

# Start behind the target in its orbit, so the radial component is negative. This constraint must always be true.
# Right as viewed from the chaser looking towards the target where orbit normal is the up direction (i.e. right = positive e_r direction in e_r/e_n plane)
function h_dock_right(t,x)
    return x[1] - Delta;
end

function hdot_dock_right(t,x)
    return x[3];
end

# Two norm velocity constraint is a single constraint
# function H_velocity(t,x)
#     H = sqrt(x[3]^2 + x[4]^2) - v_max;
# end

# Inf norm velocity constraint is effectively four constraints
function H_velocity(t,x)
    H = [x[3]; -x[3]; x[4]; -x[4]] .- v_max;
end

###################################################################################################
###################################################################################################

function f_land(t,x)
    r = x[1:3];
    v = x[4:6];
    return [v; gravity(t, r)];
end

function g_land(t,x)
    return [0.0 0.0 0.0;
            0.0 0.0 0.0;
            0.0 0.0 0.0;
            1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0];
end

function f_dock(t,x)
    return gravity(t, x);
end

function g_dock(t,x)
    return [0.0 0.0;
            0.0 0.0;
            1.0 0.0;
            0.0 1.0];
end


###################################################################################################
###################################################################################################

function H_two_norm(t,x,h,hdot,a_max;delta=0)
    function func(t,x)
        return h(t,x) + absSq(hdot(t,x) + w_xmax)/(2*a_max) - delta/a_max;
    end
    H = func(t,x);
    partial_H = ForwardDiff.derivative(y -> func(y,x), t);
    grad_H = ForwardDiff.gradient(y -> func(t,y), x)';
    global last_H = H;
    return H, partial_H, grad_H;
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

function H_docking(t,x,delta)
    H1, partial_H1, grad_H1 = H_two_norm(t,x,h_dock,hdot_dock,a_max_dock,delta=delta);
    H0l, partial_H0l, grad_H0l = H_two_norm(t,x,h_dock_left,hdot_dock_left,a_max_lr);
    H0r, partial_H0r, grad_H0r = H_two_norm(t,x,h_dock_right,hdot_dock_right,a_max_lr);

    Hv = H_velocity(t,x);
    partial_Hv = ForwardDiff.derivative(y -> H_velocity(y,x), t);
    # grad_Hv = ForwardDiff.gradient(y -> H_velocity(t,y), x)';
    grad_Hv = ForwardDiff.jacobian(y -> H_velocity(t,y), x);

    H = [H1; Hv; H0r; H0l];
    global last_H = H;
    partial_H = [partial_H1; partial_Hv; partial_H0r; partial_H0l];
    grad_H = [grad_H1; grad_Hv; grad_H0r; grad_H0l];
    return H, partial_H, grad_H;
end


###################################################################################################
###################################################################################################

function set_mu_phi(x); global mu = x; end
function set_rho(x); global rho = x; end
function set_Delta(x); global Delta = x; end
function set_u_max(x); global u_max = x; end
function set_a_max_dock(x); global a_max_dock = x; end
function set_a_max_lr(x); global a_max_lr = x; end
function set_h(x); global h = x; end
function set_hdot(x); global hdot = x; end
function set_w_umax(x); global w_umax = x;end
function set_w_xmax(x); global w_xmax = x; end
function set_gravity(x); global gravity = x; end
function set_partial_gravity(x); global partial_gravity = x; end
function set_grad_gravity(x); global grad_gravity = x; end
function set_obstacle(x); global obstacle = x; end
function set_v_max(x); global v_max = x; end

end

nothing