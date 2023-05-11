function [h,A,b] = CBF_obs(t,x,T,index,setting)
% Function for h_1 through h_4
global A_sys psi_h_nh
persistent e_max hddot_max
if isempty(e_max)
    e_max = max(svd(A_sys(3:4,:)));
end

r = x(1:2);
v = x(3:4);

[r0,v0] = obstacle_location(index,t);
rho = 1e3;
kappa = rho - norm(r - r0);
kappa_dot = -(r-r0)'*(v-v0)/norm(r-r0);

gamma = 1;
h = kappa + gamma*kappa_dot;

if nargout > 1 || (nargin==5 && isequal(setting,'nonlin'))
    if nargin==4 || isequal(setting, 'lin')
        % One can try running the simulation with this alternate bound, but it is really
        % overconservative, so we use the improved equation below instead.
        % hddot_max = e_max*local_bound(x);
        if index==1
            % This value is identical for all four obstacles, so only compute it once
            hddot_max = local_bound_improved(A_sys(3:4,:),x);
        end
    
        Lfh = -(r-r0)'/norm(r-r0);
        kappa_dot_unc = Lfh*(v-v0);

        b = -(kappa + hddot_max*(T^2/2 + gamma*T) + kappa_dot_unc*(T+gamma));
        A = Lfh*(T+gamma);
    elseif isequal(setting, 'nonlin')
        delta = T/(psi_h_nh+1);
        tspan = linspace(t, t+T, psi_h_nh+1);

        A = A_sys;
        [~, x_sim] = ode45(@(t,x) A*x, tspan, x);
        h = zeros(psi_h_nh,1);
        for i=1:psi_h_nh
            r = x_sim(i,1:2)';
            v = x_sim(i,3:4)';
            [r0,v0] = obstacle_location(index,t);
            kappa = rho - norm(r-r0);
            kappa_dot = -(r-r0)'*(v-v0)/norm(r-r0);
            % Note that for this formula to be accurate, hddot_max must be computed by
            % calling the function with the linear setting every so often
            h(i) = kappa + (gamma + delta)*kappa_dot + hddot_max*(delta^2/2 + gamma*delta);
        end
    else
        error('Unrecognized setting in CBF_obs');
    end
end
end

function bd = local_bound_improved(rows,x)
% This function is an improvement on local_bound(x) that is specific to the form of CBF
% used in this function file, and that decreases conservatism compared to
% e_max*local_bound(x). This function can still be drastically improved by considering the
% second term of hddot = -r'*u/norm(r) - norm(cross(r,v))^2/norm(r)^3. However,
% determining the maximum value of the full hddot expression is a nonconvex optimization
% that should only be considered offline.
global P
level = x'*P*x;
[v,r] = eig(P);
A0 = v';
b0 = sqrt(level)./diag(r);
A = [A0; -A0; 0, 1, 0, 0]; % docking constraint says x2 can only ever be negative
b = [b0; b0; 0];
[~,val1a] = linprog(-rows(1,:),A,b,[],[],[],[],optimset('Display','off'));
[~,val1b] = linprog(rows(1,:),A,b,[],[],[],[],optimset('Display','off'));
[~,val2] = linprog(-rows(2,:),A,b,[],[],[],[],optimset('Display','off'));
% Use mininum here, because a negative r'*u/norm(r) outcome causes h to increase
val1 = min(val1a, val1b);
bd = sqrt(val1^2 + val2^2);
end