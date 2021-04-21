%% Get Lipschitz Constants of f and g
global constants
constants.u_max = [5; 0.25];
constants.u_min = [0; -0.25];
constants.sigma = 1;
constants.rho = 10;
T = 0.1;

syms t real
g = [cos(t), 0; sin(t), 0; 0, 1];
l_g = svd(g)

%% Get l_h
func = @(x) -norm(GradH_func(x));
guess = [-10; -10; pi/4];
N = 100;
clear res fval;
for i=1:N
    guess = randn(3,1);
    try
    [res(:,i), fval(i)] = fmincon(func, guess, [], [], [], [], [], [], @constraint);
    end
end
[val, i] = min(fval);
r = res(:,i);
val
constraint(res)

l_h = 1.082;

%%
[x, y] = meshgrid(-20:1:20, -20:1:20);
phi = pi/4*ones(size(x));

dh = zeros(size(x));
for i=1:size(x,1)
    for j=1:size(x,2)
        dh(i,j) = func([x(i,j); y(i,j); phi(i,j)]);
    end
end

if 0
figure(1); clf;
surf(x, y, dh, 'EdgeAlpha', 0);
xlabel 'x'; ylabel 'y';
end

%% Get k_M
v_max = 5;
omega_max = 0.25;

k_M = three_digits(norm([5 0.25]));

%% Get max(\kappa)
func = @(x) -Hddot_func(x(1:3), x(4:5));
lower = -[inf; inf; inf; 0;     omega_max];
upper =  [inf; inf; inf; v_max; omega_max];

clear res fval;
guess = [-10; -10; pi/4; v_max; 0];
N = 100;
for i=1:N
    guess = randn(5,1);
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

kappa_max = 2.639;

%% Get global nu_1
func = @(x) -norm(GradLgH_func(x(1:3)));
clear res fval;
N = 100;
for i=1:N
    guess = randn(3,1);
    try
    [res(:,i), fval(:,i)] = fmincon(func, guess, [], [], [], [], [-40; -30; -pi], [25; 25; pi], @constraint);
    if sum(constraint(res(:,i)) > 0)
        fval(:,i) = 1;
    end
    end
end
[val, i] = min(fval);
r = res(:,i);
val
constraint(r)

l_LgH = 227.3;

nu1 = (0 + l_LgH*norm([5 0.25]) + l_h)*k_M*T
delta_nu1_min = l_h*k_M*T % increase Gamma as high as possible

l1 = 0 + l_LgH*norm([5 0.25]) + l_h;
l2 = l_LgH*norm([5 0.25]);
nu0 = l1*k_M/l2*(exp(l2*T)-1);

delta_nu0_inf = l_h*k_M/l2*(exp(l2*T)-1)

tau = [0.1, 0.01, 0.001];
delta_nu1_min_all = l_h*k_M*tau
delta_nu0_inf_all = l_h*k_M/l2*(exp(l2*tau)-1)

%%
[x, y] = meshgrid(-20:1:20, -20:1:20);
phi = pi/4*ones(size(x));

func = @(x) -norm(LgH_func(x(1:3)));
d = zeros(size(x));
for i=1:size(x,1)
    for j=1:size(x,2)
        d(i,j) = func([x(i,j); y(i,j); phi(i,j)]);
    end
end

if 0
figure(1); clf;
surf(x, y, d, 'EdgeAlpha', 0);
xlabel 'x'; ylabel 'y';
end

%% Get Global nu_2
Gamma = 1;
alpha_func = @(x) -Gamma*H_func(x);
m1_func = @(x, x0) -(LfH_func(x) - LfH_func(x) ...
    + dot(LgH_func(x) - LgH_func(x0), x(4:5)) ...
    - (alpha_func(x) - alpha_func(x0)));
func = @(x) m1_func(x(1:5), x(6:8));
clear res fval;
N = 1;
lower = [-40; -30; -pi; [0; -0.25]; -40; -30; -pi];
upper = [25; 25; pi; [5; 0.25]; 25; 25; pi];
for i=1:N
    guess = randn(8,1)*10;
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

nu2 = 0.6908;

%% Get Global \delta_nu_2 minimum
% The following optimization is nonconvex, and very dependent on the initial conditions.
% I recommend looking more closely at the plotted version.
if 0
m1_func = @(x, x0, Gamma) -(LfH_func(x) - LfH_func(x0) ...
    + dot(LgH_func(x) - LgH_func(x0), x(4:5)) ...
    - Gamma*(alpha_func(x) - alpha_func(x0)));
% Want to find the maximum of nu_2, where nu_2 = -m1_func, so we optimize
% over +m1_func.
func = @(x, Gamma) m1_func(x(1:5), x(6:8), Gamma);
clear res fval;
N = 1;
lower = [-40; -30; -pi; [0; -0.25]; -40; -30; -pi; 0];
upper = [25; 25; pi; [5; 0.25]; 25; 25; pi; Inf];
% The margin is then nu_2(Gamma)/Gamma, which we seek to minimize. The
% result of func is -nu_2, so we minimize -func/Gamma
tot_func = @(Gamma) -nested_func(Gamma, func, guess(1:8), lower(1:8), upper(1:8))/Gamma;
for i=1:N
    guess = randn(9,1)*10;
    tot_func = @(Gamma) -nested_func(Gamma, func, guess(1:8), lower(1:8), upper(1:8))/Gamma;
%     guess(9) = 100;
    try
    [res(:,i), fval(:,i)] = fmincon(tot_func, guess(9), [], [], [], [], lower(9), upper(9));
    if sum(constraint2(res(:,i)) > 0)
        fval(:,i) = 1;
    end
    end
end
[val, i] = min(fval);
[~,r0] = nested_func(val, func, guess(1:8), lower(1:8), upper(1:8));
r = [r0; res(:,i)]
val
constraint2(r)

%%
res = [];
fval = [];
slopes = [1:20:200, 2000, 4000];
for i=slopes
    m1_func = @(x, x0) -(LfH_func(x) - LfH_func(x0) ...
        + dot(LgH_func(x) - LgH_func(x0), x(4:5)) ...
        - i*(alpha_func(x) - alpha_func(x0)));
    func = @(x) m1_func(x(1:5), x(6:8));
    for j=1:10
        try
            guess = randn(8,1)*10;
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
figure(10); clf;
plot(slopes(1:length(fval)), fval./slopes(1:length(fval)));

delta_nu2_min = 0.53;
end

return
%% Set up constants
constants.lf = 0;
constants.lg = 1;
constants.n = 3;
constants.m = 2;
constants.kappa = kappa_max;
constants.nu_min = nu2;
constants.lh = l_h;
constants.kM = k_M;
constants.l_LfH_global = 0;
constants.l_LgH_global = l_LgH;

function [C, Ceq] = constraint(x)
C(1) = H_func(x);
C(2) = h0_func(x);
Ceq = [];
end

function [C, Ceq] = constraint2(x)
global constants
C(1) = H_func(x(1:3));
C(2) = h0_func(x(1:3));
C(3) = norm(x(1:2) - x(6:7)) - constants.u_max(1)*0.001;
C(4) = abs(x(3) - x(8)) - constants.u_max(2)*0.001;
C(5) = H_func(x(6:8));
C(6) = h0_func(x(6:8));
Ceq = [];
end

function [out, r] = nested_func(x, func, guess, lower, upper)
[r, out] = fmincon(@(y) func(y, x), guess, [], [], [], [], lower, upper, @constraint2);
end

function out = three_digits(x)
out = ceil(x*1e3)/1e3;
end