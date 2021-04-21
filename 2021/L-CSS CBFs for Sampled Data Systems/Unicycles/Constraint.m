function [A, b] = Constraint(x0, T, index)
% Returns the A and b constraint matrices for ensuring avoidance is achieved.
% Results should be fed into QP as A \leq b
global is_error outdata m5_margin method constants
if isempty(m5_margin)
    m5_margin = nan;
end

outdata.h(index,1) = h0_func(x0);
outdata.hdot(index,1) = 0; % no inertia in this system

H = H_func(x0);
if H > 0
    is_error = 3;
end
if outdata.h(index,1) > 0
    is_error = 2;
end
outdata.H(index,1) = H;

gamma = 1;
Gamma = 1;
persistent alpha_func m1_func m4 m1_lf m1_lg m1_lh
if isempty(alpha_func)
    alpha_func = @(x) -Gamma*H_func(x);
    
    m1_func = @(x, x0) -(LfH_func(x) - LfH_func(x) ...
        + dot(LgH_func(x) - LgH_func(x0), x(4:5)) ...
        - (alpha_func(x) - alpha_func(x0)));
    m1_lh = @(x) -Gamma*norm(GradH_func(x));
    m1_lf = @(x) -norm(GradLfH_func(x));
    m1_lg = @(x) -norm(GradLgH_func(x));
    m4 = @(x) -Hddot_func(x(1:3), x(4:5));
end
m1 = @(x) m1_func(x,x0);

Lf = LfH_func(x0);
Lg = LgH_func(x0);

u_min = constants.u_min;
u_max = constants.u_max;
[s1, s2] = reachable_set(x0, T);
    % For the unicycle system, we calculate the actual reachable set to optimize over
lower = [s1; u_min];
upper = [s2; u_max];
guess = [x0; u_max];
opts = optimset('Display', 'off');

%% Determine the margin based on all methods
% Theorem 3
[~, fval] = fmincon(m4, guess, [], [], [], [], lower, upper, [], opts);
outdata.margin_m4optim = max(-1/2*fval*T, 0);
% Corollary 3
outdata.margin_m4global = 1/2*constants.kappa*T;
% Theorem 2
[~, fval] = fmincon(m1, guess, [], [], [], [], lower, upper, [], opts);
outdata.margin_m1optim = max(-fval, 0);
% Corollary 2
outdata.margin_m1gloopt = constants.nu_min;
% Theorem 1
[~, fval1] = fmincon(m1_lh, guess, [], [], [], [], lower, upper, [], opts);
[~, fval2] = fmincon(m1_lg, guess, [], [], [], [], lower, upper, [], opts);
outdata.margin_m1local = (max(-fval1, 0) + max(-fval2, 0)*norm(u_max) + m1_lf(x0))*constants.kM*T;
% Corollary 1
outdata.margin_m1global = (constants.l_LfH_global + constants.l_LgH_global*norm(u_max) + Gamma*constants.lh)*constants.kM*T;

%% Assign Final Variables
if isequal(method, 'm1local')
    margin = outdata.margin_m1local;
    RHS = -margin - Gamma*H;
elseif isequal(method, 'm1global')
    margin = outdata.margin_m1global;
    RHS = -margin - Gamma*H;
elseif isequal(method, 'm1optim')
    margin = outdata.margin_m1optim;
    RHS = -margin - Gamma*H;
elseif isequal(method, 'm1gloopt')
    margin = outdata.margin_m1gloopt;
    RHS = -margin - Gamma*H;
elseif isequal(method, 'm4optim')
    margin = outdata.margin_m4optim;
    RHS = -margin - gamma*H/T;
elseif isequal(method, 'm4global')
    margin = outdata.margin_m4global;
    RHS = -margin - gamma*H/T;
else
    error('Unknown method type');
end

outdata.margin(index,1) = margin;
outdata.RHS(index,1) = RHS;

A = Lg;
b = RHS - Lf;

end

function [lower, upper] = reachable_set(x, T)
global constants
phi = mod(x(3), 2*pi); % it is essential that phi is positive in this method
v = constants.u_max(1);
omega = constants.u_max(2);

c1 = cos(phi - omega*T);
c2 = cos(phi + omega*T);
s1 = sin(phi - omega*T);
s2 = sin(phi + omega*T);
if sign(s1) ~= sign(s2)
    % phi is near 0 or pi, so max/minimizing x does not mean choosing the largest w
    if phi < pi/2
        % negative omega maximizes x
        omega_x1 = (0 - phi)/T;
        omega_x2 = omega;
    elseif phi < pi
        % positive omega minimizes x
        omega_x1 = -omega;
        omega_x2 = (pi - phi)/T;
    elseif phi < 3*pi/2
        % negative omega minimizes x
        omega_x1 = (pi - phi)/T;
        omega_x2 = omega;
    else
        % positive omega maximizes x
        omega_x1 = -omega;
        omega_x2 = (2*pi - phi)/T;
    end
else
    omega_x1 = -omega;
    omega_x2 =  omega;
end

if sign(c1) ~= sign(c2)
    % phi is near pi/2 or 3*pi/2, so max/minimizing y does not mean choosing the largest w
    if phi < pi/2
        % positive omega maximizes y
        omega_y1 = -omega;
        omega_y2 = (pi/2 - phi)/T;
    elseif phi < pi
        % negative omega mazimizes y
        omega_y1 = (pi/2 - phi)/T;
        omega_y2 = omega;
    elseif phi < 3*pi/2
        % positive omega minimizes y
        omega_y1 = -omega;
        omega_y2 = (3*pi/2 - phi)/T;
    else
        % negative omega minimizes y
        omega_y1 = (3*pi/2 - phi)/T;
        omega_y2 = omega;
    end
else
    omega_y1 = -omega;
    omega_y2 =  omega;
end
% There is likely redundancy in the above trees, but it is meant to be readable.

dx1 = v/omega_x1*(sin(phi + omega_x1*T) - sin(phi));
dx2 = v/omega_x2*(sin(phi + omega_x2*T) - sin(phi));
dy1 = v/omega_y1*(cos(phi) - cos(phi + omega_y1*T));
dy2 = v/omega_y2*(cos(phi) - cos(phi + omega_y2*T));

lower = [min([x(1), x(1)+dx1, x(1)+dx2]); min([x(2), x(2)+dy1, x(2)+dy2]); phi-omega*T];
upper = [max([x(1), x(1)+dx1, x(1)+dx2]); max([x(2), x(2)+dy1, x(2)+dy2]); phi+omega*T];
end