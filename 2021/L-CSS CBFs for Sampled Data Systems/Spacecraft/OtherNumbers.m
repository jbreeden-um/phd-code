% This file computes all the needed constants for these simulations. Note that many of
% these are computed as the solution to nonconvex optimizations, so the values actually
% used are given explicitly.

%% Get Lipschitz Constants of f and g
lg = 0; % analytic

global constants
u_max = 0.01; w_max = 0.2;
constants.u_max = [u_max; u_max; u_max];
constants.w_max = [w_max; w_max; w_max];
constants.mu = 0.005;
constants.theta = pi/5;
constants.s = [0;0;1];

func = @(x) -norm(Gradf_func(x));
guess = randn(6,1);
lower = -[ones(3,1); constants.w_max];
upper =  [ones(3,1); constants.w_max];
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, @constraint)

lf = 1.059;

%% Get l_h
func = @(x) -norm(GradH_func(x));
guess = randn(6,1);
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, @constraint)

lh = 58.80;

%% Get k_M
lower = -[ones(3,1); constants.w_max; constants.u_max];
upper =  [ones(3,1); constants.w_max; constants.u_max];
func = @(x) -norm(f_func(x(1:6)) + g_func(x(1:6))*x(7:9));

guess = randn(9,1);
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, @constraint)

kM = 0.347; 
% This is producing somewhat variable results, but you only need to run it
% a few times to get this.

%% Get max(\kappa)
func = @(x) -Hddot_func(x(1:6), x(7:9));
lower = -[ones(3,1); constants.w_max; constants.u_max];
upper =  [ones(3,1); constants.w_max; constants.u_max];

N = 1;
clear res fval;
for i=1:N
    guess = randn(9,1);
    try
    [res(:,i), fval(:,i)] = fmincon(func, guess, [], [], [], [], lower, upper, @constraint);
    if sum(constraint(res(:,i)) > 0)
        fval(:,i) = 1;
    end
    end
end
[val, i] = min(fval);
r = res(:,i);
val
constraint(r)

kappa_max = 2.54;
% This is producing different results depending on the initial conditions,
% which is weird since the real-time version was fairly stable. I ran it 1000 times
% and rounded up to be sure.

%% Get global nu_1
func = @(x) -norm(GradLfH_func(x(1:6)));
lower = -[ones(3,1); constants.w_max];
upper =  [ones(3,1); constants.w_max];

guess = randn(6,1);
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, @constraint)
constraint(res)

l_LfH = 22.31;

func = @(x) -norm(GradLgH_func(x(1:6)));
guess = randn(6,1);
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, @constraint)
constraint(res)

l_LgH = 217.3;

% The results were moderately consistent for both of these. You will
% probably get the correct result within 5 tries.

%% Get Global nu_2
Gamma = 1;
alpha_func = @(x) -Gamma*H_func(x);
m1_func = @(x, x0) -(LfH_func(x) - LfH_func(x) ...
    + dot(LgH_func(x) - LgH_func(x0), x(7:9)) ...
    - (alpha_func(x) - alpha_func(x0)));
func = @(x) m1_func(x(1:9), x(10:15));
clear res fval;
N = 1;
lower = -[ones(3,1); constants.w_max; constants.u_max; ones(3,1); constants.w_max];
upper =  [ones(3,1); constants.w_max; constants.u_max; ones(3,1); constants.w_max];
for i=1:N
    guess = randn(15,1);
    try
    [res(:,i), fval(:,i)] = fmincon(func, guess, [], [], [], [], lower, upper, @constraint2);
    if sum(constraint2(res(:,i)) > 0)
        fval(:,i) = 1;
    end
    end
end
[val, i] = min(fval);
r = res(:,i);
val
constraint2(r)

nu2 = 0.8815;
% Results were sort of inconsistent, but not incredibly so. I ran it for 1000 to be sure.

%% Get global delta_nu_2
% The following section computes \delta_2^g. This is a 15-dimensional optimization,
% and is very sensitive to the initial conditions. I recommend looking at the plotted
% results rather than just taking the result of the first optimization as truth.

% Also, this takes a while to run, so you may want to turn it off. Luckily, this is done
% offline, and is not even needed for simulation, so it does not matter how long it takes.
if 0
Gamma = 1;
alpha_func = @(x) -Gamma*H_func(x);
m1_func = @(x, x0) -(LfH_func(x) - LfH_func(x0) ...
    + dot(LgH_func(x) - LgH_func(x0), x(7:9)) ...
    - (alpha_func(x) - alpha_func(x0)));
func = @(x) m1_func(x(1:9), x(10:15));
clear res fval;
N = 100;
lower = -[ones(3,1); constants.w_max; constants.u_max; ones(3,1); constants.w_max];
upper =  [ones(3,1); constants.w_max; constants.u_max; ones(3,1); constants.w_max];
for i=1:N
    guess = randn(15,1);
    try
    [res(:,i), fval(:,i)] = fmincon(func, guess, [], [], [], [], lower, upper, @constraint2);
    if sum(constraint2(res(:,i)) > 0)
        fval(:,i) = 1;
    end
    end
end
[val, i] = min(fval);
r = res(:,i);
val
constraint2(r)

%%
res = [];
fval = [];
slopes = 1:20:200;
for i=slopes
    m1_func = @(x, x0) -(LfH_func(x) - LfH_func(x0) ...
        + dot(LgH_func(x) - LgH_func(x0), x(7:9)) ...
        - i*(alpha_func(x) - alpha_func(x0)));
    func = @(x) m1_func(x(1:9), x(10:15));
    for j=1:10
        try
            guess = randn(15,1);
%             guess(10:12) = -guess(1:3);
            tres = [];
            tfval = [];
            [tres(:,end+1), tfval(:,end+1)] = fmincon(func, guess, [], [], [], [], lower, upper, @constraint2);
            if sum(constraint2(tres(:,end)) > 0)
                tfval(:,end) = 1;
            end
        end
    end
    try
    [fval(:,end+1), j] = min(tfval);
    res(:,end+1) = tres(:,j);
    catch
    fval(:,end+1) = 0;
    res(:,end+1) = zeros(8,1);
    end
end
figure(7); clf;
plot(slopes(1:length(fval)), fval./slopes(1:length(fval)));

delta_nu2_min = 0.81;
end

%% Values
Gamma = 1; T = 0.1;
nu1 = (l_LfH + l_LgH*norm(constants.u_max) + lh*Gamma)*kM*T;
nu0 = (l_LfH + l_LgH*norm(constants.u_max) + lh*Gamma)*kM ...
    * (exp((l_LfH + l_LgH*norm(constants.u_max))*T) - 1)/(l_LfH + l_LgH*norm(constants.u_max));
cM = lf + lg*norm(u_max);
eM = kM*((exp(cM*T)-1)/cM - T);

m4_margin_Hdot = kappa_max*T/2
m1_margin_Hdot = nu1

m4_margin = kappa_max*T^2/2
m1_margin = nu1/Gamma

delta_nu1_inf = lh*kM*T
delta_nu0_inf = lh*kM*(exp((l_LfH + l_LgH*norm(constants.u_max))*T) - 1)/(l_LfH + l_LgH*norm(constants.u_max))

tau = [0.1, 0.01, 0.001];
delta_nu1_inf_all = lh*kM*tau
delta_nu0_inf_all = lh*kM*(exp((l_LfH + l_LgH*norm(constants.u_max))*tau) - 1)/(l_LfH + l_LgH*norm(constants.u_max))

return

%% Setting values
% Run this section once. I've put it after the return statement so you can play around
% with this file without breaking the simulations.
constants.lf = lf;
constants.lg = lg;
constants.n = 6;
constants.m = 3;
constants.kappa = kappa_max;
constants.lh = lh;
constants.kM = kM;
constants.l_LfH_global = l_LfH;
constants.l_LgH_global = l_LgH;
constants.nu_min = nu2;
constants.nu0 = nu0;

function [C, Ceq] = constraint(x)
C(1) = H_func(x);
C(2) = h0_func(x);
Ceq = norm(x(1:3)) - 1;
end

function [C, Ceq] = constraint2(x)
global constants
C(1,1) = H_func(x(1:6));
C(2,1) = h0_func(x(1:6));
C(3,1) = H_func(x(10:15));
C(4,1) = h0_func(x(10:15));
C(5:7,1) = abs(x(1:3) - x(10:12)) - constants.w_max*0.001 - constants.u_max.*0.001^2;
    % I deliberately did not divide the last term by 2 because attitude
    % dynamics are weird and I want there to be margin
C(8:10,1) = abs(x(4:6) - x(13:15)) - constants.u_max*0.001;
Ceq = [];
end