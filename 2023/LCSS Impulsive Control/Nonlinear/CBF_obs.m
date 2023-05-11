function [h,A,b] = CBF_obs(t,x,T,index,setting)
% Function for h_1 through h_4
global psi_h_nh robust_margin_checking
persistent hddot_max

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
        hddot_max = 7.59e-7*sqrt(compute_V(t,x));
    
        Lfh = -(r-r0)'/norm(r-r0);
        kappa_dot_unc = Lfh*(v-v0);

        b = -(kappa + hddot_max*(T^2/2 + gamma*T) + kappa_dot_unc*(T+gamma));
        A = Lfh*(T+gamma);
    elseif isequal(setting, 'nonlin')
        delta = T/psi_h_nh;
        tspan = linspace(t, t+T, psi_h_nh+1);

        [~, x_sim] = ode45(@two_body, tspan, x);
        h = zeros(psi_h_nh,1);
        
        if robust_margin_checking
            V = 0;
            for i=1:length(tspan)
                V = max(V, compute_V(tspan(i), x_sim(i,:)));
            end
            hddot_max = 7.59e-7*sqrt(V);
        end
        
        for i=1:psi_h_nh
            r = x_sim(i,1:2)';
            v = x_sim(i,3:4)';
            [r0,v0] = obstacle_location(index,tspan(i));
            kappa = rho - norm(r-r0);
            kappa_dot = -(r-r0)'*(v-v0)/norm(r-r0);
            % Note that for this formula to be accurate, hddot_max must be computed by
            % calling the function with the linear or robust setting every so often
            h(i) = kappa + (gamma + delta)*kappa_dot + hddot_max*(delta^2/2 + gamma*delta);
        end
    else
        error('Unrecognized setting in CBF_obs');
    end
end
end