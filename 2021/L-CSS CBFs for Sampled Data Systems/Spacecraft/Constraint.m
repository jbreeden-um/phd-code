function [A, b] = Constraint(x0, T, index)
% Returns the A and b constraint matrices for ensuring avoidance is achieved.
% Results should be fed into QP as A \leq b
global is_error outdata constants method

outdata.h(index,1) = h0_func(x0);
outdata.hdot(index,1) = h0dot_func(x0);

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
        + dot(LgH_func(x) - LgH_func(x0), x(7:9)) ...
        - (alpha_func(x) - alpha_func(x0)));
    m1_lh = @(x) -Gamma*norm(GradH_func(x));
    m1_lf = @(x) -norm(GradLfH_func(x));
    m1_lg = @(x) -norm(GradLgH_func(x));
    
    m4 = @(x) -Hddot_func(x(1:6), x(7:9));
end
m1 = @(x) m1_func(x,x0);

Lf = LfH_func(x0);
Lg = LgH_func(x0);
u_max = constants.u_max;

delta = DeltaMax(x0, T);
lower = [x0 - delta*T; -u_max];
upper = [x0 + delta*T;  u_max];
opts = optimset('Display', 'off');
nonlin = @(x) dist_constraint(x(1:6), x0, delta*T);
    % Using the ball constraint for the spacecraft case rather than attempting to
    % figure out what the reachable set really is.

ug = u_max - dot(u_max, x0(1:3))/dot(x0(1:3), x0(1:3))*x0(1:3);
if cross(ug, x0(1:3)) > 0
    guess = [x0; ug];
else
    guess = [x0; -ug];
end
    % The following optimizations are not necessarily convex. With the above guesses, they
    % converge to the correct solution >99% of time steps.

%% Determine the margin based on all methods
% All the methods run in parallel here, so we can do comparisons, but if you are not
% interested in that, I highly recommend commenting out the unneeded methods.

% Determine the solution according to theorem 3
[~, fval] = fmincon(m4, guess, [], [], [], [], lower, upper, nonlin, opts);
outdata.margin_m4optim = max(-1/2*fval*T, 0);
% Determine the solution according to corollary 3
outdata.margin_m4global = 1/2*constants.kappa*T;

% Determine the solution according to theorem 2
[~, fval] = fmincon(m1, guess, [], [], [], [], lower, upper, nonlin, opts);
outdata.margin_m1optim = max(-fval, 0);
% Determine the solution according to corollary 2
outdata.margin_m1gloopt = constants.nu_min;

% Determine the solution according to theorem 1
[~, fval1] = fmincon(m1_lf, guess, [], [], [], [], lower, upper, [], opts);
[~, fval2] = fmincon(m1_lg, guess, [], [], [], [], lower, upper, [], opts);
[~, fval3] = fmincon(m1_lh, guess, [], [], [], [], lower, upper, [], opts);
outdata.margin_m1local = (max(-fval1, 0) + max(-fval2, 0)*norm(u_max) + max(-fval3, 0))*delta*T;
% Determine the solution according to corollary 1
outdata.margin_m1global = (constants.l_LfH_global + constants.l_LgH_global*norm(u_max) + Gamma*constants.lh)*constants.kM*T;
% Determine the solution according to lemma 2
outdata.margin_m0 = constants.nu0;

%% Assign Final Variables
if isequal(method, 'm0')
    margin = outdata.margin_m0;
    RHS = -margin - Gamma*H;
elseif isequal(method, 'm1local')
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

function out = DeltaMax(x, T)
global constants
u_max = constants.u_max;
lf = constants.lf;
lg = constants.lg;
num = norm(f_func(x)) + norm(g_func(x))*norm(u_max);
den = 1 - (lf + lg*norm(u_max))*T;
if den <= 0
    disp 'Invalid Approximation. Choose a smaller time step';
end
out = num/den;
end

function [C, Ceq] = dist_constraint(x, x0, delta)
global constants
n = constants.n; % dimension of the system
num = round(length(x)/n);
C = zeros(num,1);
Ceq = zeros(num,1);
for i=1:num
    C(i) = norm(x((n*(i-1)+1):(n*i)) - x0((n*(i-1)+1):(n*i))) - delta;
    Ceq(i) = norm(x((6*(i-1)+1):(6*(i-1)+3))) - 1;
end
end
