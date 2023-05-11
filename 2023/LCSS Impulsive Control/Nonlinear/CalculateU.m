function [u, compute] = CalculateU(t,x,T,dt)
global use_dock_constraint
if sum(isnan(x))
    disp('x is NaN');
end
    
global psi_v_min_rate psi_v_rhs_max psi_v_tol_dt
V = compute_V(t,x);
c = zeros(6,1,'logical');
if psi_v(t+dt,t,x,psi_v_tol_dt) <= min(psi_v_rhs_max, psi_v_min_rate*V)
    c(1) = 1;
end
% Always use the linear bounds for small dt
for i=1:4
    [~, ~, b0] = CBF_obs(t,x,dt,i,'lin');
    if b0 >= 0 % if applying zero control is safe
        c(1+i) = 1;
    end
end
if use_dock_constraint
    [~,~,b0] = CBF_dock(t,x,dt,'lin');
    if b0 >= 0 % if applying zero control is safe
        c(6) = 1;
    end
else
    c(6) = 1;
end
if min(c) == 1 % if all constraints above were satisfied for zero control input
    u = [0;0];
    compute = NaN;
else
    [u, compute] = u_star(t,x,T);
end
end

function [u, compute] = u_star(t,x,T)
% Constraints on the rate of decrease of V
global psi_v_tol_dT gamma_V
tol = min(0.01*compute_V(t,x), psi_v_tol_dT);
c1 = @(u) psi_v_constraint(u, t, x, T, compute_V(t,x), tol);
c2 = @(u) V_constraint(u, x, t, gamma_V);
% The following code mimics the final cases of the main theorem, but is not necessary here
% c2 = @(u) V_explicit_constraint(u, x, t, gamma, T);

% Constraints on safety
global use_dock_constraint use_psi_h_star
A_cbf = zeros(5,2);
b_cbf = zeros(5,1);
for i=1:4
    [~,A_cbf(i,:),b_cbf(i)] = CBF_obs(t,x,T,i,'lin');
end
if use_dock_constraint
    [~,A_cbf(5,:),b_cbf(5)] = CBF_dock(t,x,T,'lin');
end
if ~use_psi_h_star
    A = A_cbf;
    b = b_cbf;
	c3 = @(u) [];
else
    A = zeros(5,2);
    b = ones(5,1);
    nlcon1 = @(u) CBF_obs(t,UpdateX_Jump(x,u(1:2)),T,1,'nonlin');
    nlcon2 = @(u) CBF_obs(t,UpdateX_Jump(x,u(1:2)),T,2,'nonlin');
    nlcon3 = @(u) CBF_obs(t,UpdateX_Jump(x,u(1:2)),T,3,'nonlin');
    nlcon4 = @(u) CBF_obs(t,UpdateX_Jump(x,u(1:2)),T,4,'nonlin');
    if use_dock_constraint
        nlcon5 = @(u) CBF_dock(t,UpdateX_Jump(x,u(1:2)),T,'nonlin');
        c3 = @(u) stack_constraints_int(u, nlcon1, nlcon2, nlcon3, nlcon4, nlcon5);
    else
        c3 = @(u) stack_constraints_int(u, nlcon1, nlcon2, nlcon3, nlcon4);
    end
end
nlcon = @(u) stack_constraints(u, c1, c2, c3);
J = eye(2);
J_slack = 10;
u_lin = quadprog(J,[0;0],A_cbf,b_cbf,[],[],[],[],[],optimset('Display','off'));
u0 = [u_lin; 0];
tic
[u_all, ~, flag] = fmincon(@(u) u(1:2)'*J*u(1:2)/2 + J_slack*u(3)^2, u0, [A, zeros(5,1)], b, ...
    [], [], [], [], nlcon, optimset('Display', 'off', 'Algorithm', 'sqp'));
compute = toc;
if flag<1
    if flag == 0
        disp(['fmincon took too many iterations, but found a solution, code = ' num2str(flag) ', t = ' num2str(t)]);
    else
        disp(['fmincon failed to find a feasible point, code = ' num2str(flag) ', t = ' num2str(t)]);
    end
end
check = max(nlcon(u_all));
if check > 0.1
    disp(['Computation error at t = ' num2str(t), ', tol = ' num2str(check)])
end
if sum(isnan(u_all))
    disp('u is NaN');
end
u = u_all(1:2);
end

function C = psi_v_constraint(u,t,x,T,V,tol)
global psi_v_min_rate psi_v_rhs_max
x_new = UpdateX_Jump(x,u(1:2));
d = u(3);
C = psi_v(t+T,t,x_new,tol,x) - min(psi_v_rhs_max, psi_v_min_rate*V) + d;
end

function C = V_constraint(u,x,t,gamma)
x_new = UpdateX_Jump(x,u(1:2));
d = u(3);
C = compute_V(t,x_new) - (1-gamma)*compute_V(t,x) + d;
end

% function C = V_explicit_constraint(u,x,t,gamma,T)
% global A_sys
% x_new = UpdateX_Jump(x,u);
% C = compute_V(t+T,expm(A_sys*T)*x_new) - (1-gamma)*compute_V(t,x);
% end

function [C, Ceq] = stack_constraints(u, c1, c2, c3)
C = [c1(u); c2(u); c3(u)];
Ceq = [];
end

function C = stack_constraints_int(u, c1, c2, c3, c4, c5)
if nargin==5
    C = [c1(u); c2(u); c3(u); c4(u)];
elseif nargin==6
    C = [c1(u); c2(u); c3(u); c4(u); c5(u)];
else
    error('Unexpected number of CBFs');
end
end