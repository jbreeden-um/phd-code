function [u, compute] = CalculateU(t,x,rhohat,Ta,Ts,Tm,sigma_m)
if sum(isnan(x))
    disp('x is NaN');
end
    
global psi_v_min_rate psi_v_rhs_max psi_v_tol_dt
V = compute_V(t,x);
c = zeros(8,1,'logical');
c(1) = psi_v(t+Ts,t,x,psi_v_tol_dt) <= min(psi_v_rhs_max, psi_v_min_rate*V);
% Always use the linear bounds for small dt
for i=1:6
    [~, ~, b0] = CBF_obs(t,x,Ts,rhohat,0,i,'lin');
    c(i+1) = b0 >= 0; % if applying zero control is safe
end
[~,~,b0] = CBF_dock(t,x,Ts,rhohat,0,'lin');
c(8) = b0 >= 0; % if applying zero control is safe
if min(c) == 1 && sigma_m > max(Ta,Tm)+Ts % if all constraints above were satisfied for zero control input
    u = [0;0];
    compute = NaN;
else
    [u, compute] = u_star(t,x,rhohat,Ta);
end
end

function [u, compute] = u_star(t,x,rhohat,T)
% Constraints on the rate of decrease of V (prior paper)
global psi_v_tol_dT gamma_V
tol = min(0.01*compute_V(t,x), psi_v_tol_dT);
c1 = @(u) psi_v_constraint(u, t, x, T, compute_V(t,x), tol);
c2 = @(u) V_constraint(u, x, t, gamma_V);

% Constraints on safety (updated for robustness)
A_cbf = zeros(7,2);
b_cbf = zeros(7,1);
for i=1:6
    [~,A_cbf(i,:),b_cbf(i)] = CBF_obs(t,x,T,rhohat,0,i,'lin');
    % This is just to get a better initial guess for the real control computation, so it
    % is okay to assume that the control input is zero.
end
[~,A_cbf(7,:),b_cbf(7)] = CBF_dock(t,x,T,rhohat,0,'lin');

A = zeros(5,2);
b = ones(5,1);
nlcon1 = @(u) CBF_obs(t,UpdateX_Jump(t,x,u(1:2)),T,rhohat,norm(u(1:2)),1,'nonlin');
nlcon2 = @(u) CBF_obs(t,UpdateX_Jump(t,x,u(1:2)),T,rhohat,norm(u(1:2)),2,'nonlin');
nlcon3 = @(u) CBF_obs(t,UpdateX_Jump(t,x,u(1:2)),T,rhohat,norm(u(1:2)),3,'nonlin');
nlcon4 = @(u) CBF_obs(t,UpdateX_Jump(t,x,u(1:2)),T,rhohat,norm(u(1:2)),4,'nonlin');
nlcon5 = @(u) CBF_obs(t,UpdateX_Jump(t,x,u(1:2)),T,rhohat,norm(u(1:2)),5,'nonlin');
nlcon6 = @(u) CBF_obs(t,UpdateX_Jump(t,x,u(1:2)),T,rhohat,norm(u(1:2)),6,'nonlin');
nlcon7 = @(u) CBF_dock(t,UpdateX_Jump(t,x,u(1:2)),T,rhohat,norm(u(1:2)),'nonlin');
    
nlcon = @(u) stack_constraints(u, c1, c2, nlcon1, nlcon2, nlcon3, nlcon4, nlcon5, nlcon6, nlcon7);
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
x_new = UpdateX_Jump(t,x,u(1:2));
d = u(3);
C = psi_v(t+T,t,x_new,tol,x) - min(psi_v_rhs_max, psi_v_min_rate*V) + d;
end

function C = V_constraint(u,x,t,gamma)
x_new = UpdateX_Jump(t,x,u(1:2));
d = u(3);
C = compute_V(t,x_new) - (1-gamma)*compute_V(t,x) + d;
end

function [C, Ceq] = stack_constraints(u, c1, c2, c3, c4, c5, c6, c7, c8, c9)
if nargin==5
    C = [c1(u); c2(u); c3(u); c4(u)];
elseif nargin==6
    C = [c1(u); c2(u); c3(u); c4(u); c5(u)];
elseif nargin==7
    C = [c1(u); c2(u); c3(u); c4(u); c5(u); c6(u)];
elseif nargin==8
    C = [c1(u); c2(u); c3(u); c4(u); c5(u); c6(u); c7(u)];
elseif nargin==9
    C = [c1(u); c2(u); c3(u); c4(u); c5(u); c6(u); c7(u); c8(u)];
elseif nargin==10
    C = [c1(u); c2(u); c3(u); c4(u); c5(u); c6(u); c7(u); c8(u); c9(u)];
else
    error('Unexpected number of CBFs');
end
Ceq = [];
end