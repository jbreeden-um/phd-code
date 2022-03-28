function [u, H] = CalculateU_MPC(t,x)
persistent x0
h = h_func(t,x);
H = h;

k = 1;
vdes = 12;

% Results are observed to only be stable when we use the same MPC time step as simulation time step
dt = .1;
% Horizon is the same as the horizon of Hstar
N = round(2.5/dt)+1;

R = [1; 1];
Q = [k; k];
cost_func = @(x,u) sum(sum(([0, 1, 0, 0; 0, 0, 0, 1]*reshape(x, [4, N]) - vdes).^2, 2).*Q) ...
    + sum(sum(reshape(u.^2, [2, N-1]), 2).*R);

cost = @(x) cost_func(x(1:4*N), x((4*N+1):end));
constraint = @(y) constraint_func(y(1:4*N), y((4*N+1):end), dt, N, x);

if t==0
    x0 = [reshape(x.*ones(1,N), [4*N, 1]); zeros(2*(N-1),1)];
end
lb = [-inf*ones(4*N,1); -10*ones(2*(N-1),1)];
ub = [inf*ones(4*N,1); 10*ones(2*(N-1),1)];
[sol, fval, flag] = fmincon(cost, x0, [], [], [], [], lb, ub, constraint, optimset('Display', 'off', 'MaxFunEval', 2e4));
if flag == 0
    disp(['t = ' num2str(t) ': fmincon is unhappy. Expect discontinuous results.']);
elseif flag < 0
    disp(['t = ' num2str(t) ': fmincon is very unhappy. Expect constraint violations.']);
end

x0 = [sol(5:4*N); zeros(4,1); sol(4*N+3:end); zeros(2,1)];
u = sol((4*N+1):(4*N+2));
end

function [C, Ceq] = constraint_func(x,u,dt,N,x_curr)
Ceq = zeros(4*N,1);
C = zeros(N-1,1);
for i=1:(N-1)
    Ceq((1+(i-1)*4):(i*4)) = UpdateX(0, dt, x((1+(i-1)*4):(i*4)), u((1+(i-1)*2):(i*2))) - x((1+i*4):(1+i)*4);
    C(i) = h_func(0, x((1+i*4):(1+i)*4));
end
Ceq(end-3:end) = x(1:4) - x_curr;
end