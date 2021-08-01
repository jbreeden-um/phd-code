function [A, b] = Constraint_AvoidanceInfNorm(x)
% Returns the A and b constraint matrices for ensuring avoidance is achieved.
% Results should be fed into QP as A \leq b
global is_error outdata
persistent a_max
if isempty(a_max)
    file = load('InData/SimData.mat');
    a_max = file.a_max;
end

rvec = x(1:3);
vvec = x(4:6);
r = norm(rvec);
rho = 11;
outdata.rc = rvec/norm(rvec)*10;

h = rho - r;

outdata.h = h;

if h > 0
    is_error = 1;
end

[t_prop, x_prop] = ode45(@odefunc, [0 100], [rvec; vvec], odeset('RelTol', 1e-8, 'Events', @stop_event), a_max);
h_prop = h_func_vec(x_prop, rho);
[H, index] = max(h_prop);

t_end = t_prop(index); % For this simple problem, the maximizer time is always unique, so we ignore the "max" in Equation 13.
x_end = x_prop(index,:)';

% One problem with this methodology is finding the true maximizer time. Generally, the
% above ode propagation will return points on either side of the true maximizer, but not
% the exact maximizer. Thus, letting H = max(h_prop) will always be an underestimate. We
% remedy this by taking the three points closest to the maximizer, and calculating a
% quadratic approximation of the curve, which in our experience always yields a good
% enough estimate for the true maximizer value.
if index==1, index = index+1; end
t_vals = t_prop((index-1):(index+1));
h_vals = h_prop((index-1):(index+1));
poly_abc = vander(t_vals)\h_vals; % quadratic polynomial coefficients
poly_a = poly_abc(1); poly_b = poly_abc(2);
t_end = max(0, -poly_b/(2*poly_a)); % unique critical point of a quadratic polynomial
x_end = interp1(t_prop, x_prop, t_end)'; % interpolate x as well, though this is less sensitive to the above error than H is
H = polyval(poly_abc, t_end); % calculate maximizer of h_prop is computed continuously
 
if H > 0
    is_error = 1;
end
outdata.H = H;

% theta(t_{c,0}) = [eye(3), eye(3)*t_{c,0}; zeros(3), eye(3)];
LgH = -x_end(1:3)'/norm(x_end(1:3))*t_end;
LfH = -x_end(1:3)'/norm(x_end(1:3))*vvec;

alpha = -H; % straight line of slope 1

A = LgH;
b = alpha - LfH;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = h_func_vec(x, rho)
r = x(:,1:3);
out = rho - vecnorm(r,2,2);
end

function out = hdot_func(x)
r = x(1:3);
v = x(4:6);
out = -dot(r,v)/norm(r);
end

function out = hddot_func(x, u)
r = x(1:3);
v = x(4:6);
out = -norm(cross(r,v))^2/norm(r)^3 - dot(r,u)/norm(r);
end

function xdot = odefunc(t, x, a_max)
rvec = x(1:3);
vvec = x(4:6);
u = a_max*sign(rvec);
vvecdot = u;
xdot = [vvec; vvecdot];
end

function [value, isterminal, direction] = stop_event(t, x, a_max)
rvec = x(1:3);
u = a_max*sign(rvec);
c2 = hddot_func(x, u);
if c2 > 0
    value = c2;
else
    value = hdot_func(x) + 100; % ensure there is enough to curve fit
end
isterminal = 1;
direction = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Above, the ODE for theta was a time invariant linear ODE, so we calculated the solution
% explicitly: theta = expm([zeros(3), t*eye(3); zeros(3), zeros(3)]);
% The following code is provided as an outline of how one would go about numerically
% implementing the same result. We also refer the reader to the code "Constraints.jl"
% written in Julia in the 2021/Automatica folder. The code there is written for a much
% more generic problem, and Julia code often runs 10x faster than MATALB.

function [out, H_val, z_val, t_end] = grad_H(x)
[H_val, t_end] = H_func(x);
y0 = [x; reshape(eye(6), [36 1])]; % put z and theta into a single state vector
opts = [];
if t_end==0
    y = y0'; % no propagation required in this case 
elseif t_end > 0
    [t, y] = ode45(@grad_odefunc, [0 t_end], y0, opts); % propagate both z and theta
else
    error('Invalid t. Check H_func');
end

z_val = y(end,1:6);
theta = reshape(y(end,7:42), [6 6]); % reshape the state vector into a matrix
r = z_val(1:3);
dh_dx = [-r/norm(r), 0, 0, 0]; % IMPORTANT: dh_dx must be evaluated at r(t_{c,0}), not r(t)
out = dh_dx*theta;
end

function xdot = grad_odefunc(t, x)
z = x(1:6);
theta = reshape(x(7:42), [6 6]);
zdot = odefunc(t,z);
Df = [zeros(3), eye(3); zeros(3), zeros(3)]; % gradients of f entries with x
% g(x) is a constant, so it does not depend on x, but u might depend on x
g = [zeros(3); eye(3)];
Dg_u = g*grad_u(z); % gradients of g*u entries with x
thetadot = (Df + Dg_u)*theta;
xdot = [zdot; reshape(thetadot, [36 1])];
end

function out = grad_u(x)
out = zeros(3, 6); % u = sign(r) is constant almost everywhere
end