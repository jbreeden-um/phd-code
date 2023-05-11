function PrePlanner

mu = 398600e9;
Re = 6378e3;
alt = 600e3;
n = sqrt(mu/(Re+alt)^3);
A_sys = [0,     0, 1,    0;
         0,     0, 0,    1;
         3*n^2, 0, 0,    2*n;
         0,     0, -2*n, 0];
dT = 300;
Phi = expm(A_sys*dT);
B_sys = [0, 0; 0, 0; 1, 0; 0, 1];
x0_hcw = [-7.4; -10e3; 0; -1];
xf_hcw = [0;0;0;0];
x0 = x0_hcw + get_center(0) + [0; 0; [eye(2), zeros(2,1)]*cross([0;0;n], [x0_hcw(1:2); 0])];

duration = ceil(3000/dT)*dT;
dim_x = 4;
dim_u = 2;
Ns = duration/dT;
N = dim_x*(Ns+1) + dim_u*Ns;

% Dynamics constraints
Aeq = zeros(dim_x*(Ns+2),N);
beq = zeros(dim_x*(Ns+2),1);
for i=1:Ns
    i1 = (i-1)*dim_x+1;
    i2 = i*dim_x;
    ip1 = i*dim_x+1;
    ip2 = (i+1)*dim_x;
    u1 = (Ns+1)*dim_x + (i-1)*dim_u + 1;
    u2 = (Ns+1)*dim_x + i*dim_u;
    Aeq(i1:i2,ip1:ip2) = -eye(dim_x);
    Aeq(i1:i2,i1:i2) = Phi;
    Aeq(i1:i2,u1:u2) = Phi*B_sys;
end

% Start and end constraints
i1 = Ns*dim_x+1;
i2 = (Ns+1)*dim_x;
Aeq(i1:i2,1:dim_x) = eye(dim_x);
beq(i1:i2,1) = x0_hcw;
i1 = (Ns+1)*dim_x+1;
i2 = (Ns+2)*dim_x;
i3 = Ns*dim_x+1;
i4 = (Ns+1)*dim_x;
Aeq(i1:i2,i3:i4) = eye(dim_x);
beq(i1:i2,1) = xf_hcw;

%% First Plot
figure(10); clf;
theta = linspace(0, 2*pi, 100);
cx = cos(theta);
cy = sin(theta);
rho = 1;
xc = get_center(0);
obs = obstacle_location(1,0)-xc(1:2); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r'); hold on;
obs = obstacle_location(2,0)-xc(1:2); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(3,0)-xc(1:2); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(4,0)-xc(1:2); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
xlabel 'x_1 (km)'; ylabel 'x_2 (km)';

Q = zeros(N);
iu1 = (Ns+1)*dim_x + 1;
iu2 = (Ns+1)*dim_x + dim_u*Ns;
Q(iu1:iu2,iu1:iu2) = eye(Ns*dim_u);
A = zeros(Ns+1,N);
b = zeros(Ns+1,1);
for i=1:(Ns+1), A(i,dim_x*(i-1)+2) = 1; end
sol_q = quadprog(Q, zeros(N,1), A, b, Aeq, beq);

x = reshape(sol_q(1:((Ns+1)*dim_x)), [dim_x, Ns+1]);
plot(x(1,:)/1e3, x(2,:)/1e3, 'm', 'LineWidth', 2)

figure(12); clf;
plot(dT*(0:1:(Ns-1)), sol_q((Ns+1)*dim_x+(1:dim_u:Ns*dim_u)), 'o'); hold on;
plot(dT*(0:1:(Ns-1)), sol_q((Ns+1)*dim_x+(2:dim_u:Ns*dim_u)), 'o');
xlabel('Time (s)'); ylabel('Control (m/s');


u = reshape(sol_q(iu1:iu2), [2 Ns]);
y = zeros(dim_x,Ns+1);
y(:,1) = x0;
for i=1:Ns
    [~,theta] = get_center(dT*(i-1));
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    x_new = UpdateX_Jump(y(:,i), R*u(:,i));
    [~,z] = UpdateX_Flow(0,x_new,dT,3);
    y(:,i+1) = z(:,end)';
end

q_cost = sol_q'*Q*sol_q/2
real_cost = sum(vecnorm(reshape(sol_q(iu1:iu2), [2, Ns])))
final_V = compute_V(duration,y(:,end))

figure(15); clf;
y_hcw = convert_to_hcw(0:dT:duration, y);
plot(y_hcw(1,:)/1e3, y_hcw(2,:)/1e3);
return

% Conclusion: Whether to use this or not really comes back to the
% "cluttered environment" question. The initial solution is safe, but I can
% change the obstacle locations so as to become unsafe. In that
% case, should we manually tune the optimizer to avoid the unsafe states,
% or do we turn to an explicit solution like fmincon?
% If fmincon, then the solution is going to take a long time to run, as illustrated if
% nlcon is enabled below.

%% Run the solver
iu1 = (Ns+1)*dim_x + 1;
iu2 = (Ns+1)*dim_x + dim_u*Ns;
% func = @(x) sum(vecnorm(reshape(x(iu1:iu2), [2, Ns])))/2;
func = @(x) x'*Q*x/2;
nlcon = []; % @(x) nl_constraints(x, Ns);

q_cost = func(sol_q)
init_cost = func(x_guess)
tic
[sol, val, flag] = fmincon(func, x_guess, [], [], Aeq, beq, [], [], nlcon, optimset('algorithm', 'sqp', 'MaxFunEvals', 2e5, 'MaxIter', 4e3));
compute = toc
func(sol)
flag

figure(10);
x = reshape(sol(1:(iu1-1)), [4, Ns+1]);
plot(x(1,:)/1e3, x(2,:)/1e3, 'k', 'LineWidth', 2)

figure(12);
plot(dT*(0:1:(Ns-1)), sol((Ns+1)*dim_x+(1:dim_u:Ns*dim_u)));
plot(dT*(0:1:(Ns-1)), sol((Ns+1)*dim_x+(2:dim_u:Ns*dim_u)));

% The solution to an identical problem does not even converge using fmincon.
% fmincon is then many times slower if we add in the nonlinear constraints as well.

end

% Note that the following nonlinear constraints should be used with the CBF_obs file in
% the Linear folder.
function [C, Ceq] = nl_constraints(x, Ns)
C = zeros((Ns+1)*5,1);
Ceq = [];
dim_x = 4;
for i=1:(Ns+1)
    i1 = (i-1)*5;
    ix0 = (i-1)*dim_x;
    xi = x((ix0+1):(ix0+dim_x));
    C(i1+1) = CBF_obs(0,xi,[],1);
    C(i1+2) = CBF_obs(0,xi,[],2);
    C(i1+3) = CBF_obs(0,xi,[],3);
    C(i1+4) = CBF_obs(0,xi,[],4);
    C(i1+5) = xi(2); %CBF_dock(0,xi,[]);
end
end


