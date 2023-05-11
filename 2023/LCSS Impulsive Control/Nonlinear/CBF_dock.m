function [h,A,b] = CBF_dock(t,x,T,setting)
% Function for h_5
global A_sys psi_h_nh
persistent e_max hddot_max
if isempty(e_max)
    e_max = max(svd(A_sys(3,:)));
end

r = x(1:2);
v = x(3:4);

kappa = r(2);
kappa_dot = v(2);

gamma = 1;
h = kappa + gamma*kappa_dot;

if nargout > 1 || (nargin==4 && isequal(setting,'nonlin'))
    if nargin==3 || isequal(setting, 'lin')
        % One can try running the simulation with this alternate bound, but it is really
        % overconservative, so we used the improved equation below instead.
%         hddot_max = e_max*local_bound(x);
        hddot_max = local_bound_improved(A_sys(4,:),x);

        Lfh = [0, 1];
        kappa_dot_unc = Lfh*v;

        b = -(kappa + hddot_max*(T^2/2 + gamma*T) + kappa_dot_unc*(T+gamma));
        A = Lfh*(T+gamma);
    elseif isequal(setting, 'nonlin')
        delta = T/(psi_h_nh+1);
        tspan = linspace(t, t+T, psi_h_nh+1);

        A = A_sys;
        [~, x_sim] = ode45(@(t,x) A*x, tspan, x);
        h = zeros(psi_h_nh,1);
        for i=1:psi_h_nh
            kappa = x_sim(i,2);
            kappa_dot = x_sim(i,4);
            % Note that for this formula to be accurate, hddot_max must be computed by
            % calling the function with the linear setting every so often
            h(i) = kappa + (gamma + delta)*kappa_dot + hddot_max*(delta^2/2 + gamma*delta);
        end
    else
        error('Unrecognized setting in CBF_dock');
    end

end
end

function bd = local_bound_improved(F,x)
% This function is an improvement on local_bound(x) that is specific to the form of CBF
% used in this function file, and that decreases conservatism compared to
% e_max*local_bound(x).
global P
level = x'*P*x;
[v,r] = eig(P);
A0 = v';
b0 = sqrt(level)./diag(r);
A = [A0; -A0; 0, 1, 0, 0]; % docking constraint says that x2 can only ever be negative
b = [b0; b0; 0];
[~,val] = linprog(-F,A,b,[],[],[],[],optimset('Display','off'));
bd = -val;
end