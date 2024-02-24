function [h,A,b,hhat] = CBF_obs(t,x,T,rhohat,u_norm,index,setting)
% Function for h_1 through h_6
global psi_h_nh
persistent hddot_max

r = x(1:2);
v = x(3:4);

x0 = obstacle_location(index,t);
r0 = x0(1:2);
v0 = x0(3:4);
rho = 1e3;
kappa = rho - norm(r - r0);
kappa_dot = -(r-r0)'*(v-v0)/norm(r-r0);

gamma = 1;
h = kappa + gamma*kappa_dot;

% The following Lipschitz constants could be generalized to the trajectory more precisely,
% but instead we simply pad them by 2% to account for variation over the trajectory.
if nargout == 4 && nargin==7 && isequal(setting,'post')
    l_margin = 1;
    A = [];
    b = [];
else
    l_margin = 1.02;
end
l_r = 1 + l_margin*2*gamma*norm(v-v0)/norm(r-r0);
l_v = gamma;

hhat = h + l_r*rhohat(1) + l_v*rhohat(2);

if nargout == 2 || nargout == 3 || (nargin==7 && isequal(setting,'nonlin'))
    if nargin==6 || isequal(setting, 'lin')
        hddot_max = 7.6e-7*sqrt(compute_V(t,x));
    
        Lfh = -(r-r0)'/norm(r-r0);
        kappa_dot_unc = Lfh*(v-v0);
        rho_sim = UpdateEstimates_Flow(t,x,rhohat,T);
        
        b = -(kappa + hddot_max*(T^2/2 + gamma*T) + kappa_dot_unc*(T+gamma) + l_r*rho_sim(1) + l_v*rho_sim(2));
        A = Lfh*(T+gamma);
    elseif isequal(setting, 'nonlin')
        delta = T/psi_h_nh;

        [~, rhohat] = UpdateEstimates_Jump(t,x,rhohat,[u_norm;0]);
        [tspan, x_sim, rho_sim] = UpdateEstimates_Flow(t,x,rhohat,T,psi_h_nh+1);
        h = zeros(psi_h_nh,1);
        
        for i=1:psi_h_nh
            r = x_sim(1:2,i);
            v = x_sim(3:4,i);
            rho_r = rho_sim(1,i+1);
            rho_v = rho_sim(2,i+1);
            x0 = obstacle_location(index,tspan(i));
            r0 = x0(1:2); v0 = x0(3:4);
            kappa = rho - norm(r-r0);
            kappa_dot = -(r-r0)'*(v-v0)/norm(r-r0);
            
            % Note that for this formula to be accurate, hddot_max must be computed by
            % calling the function with the linear or robust setting every so often.
            %
            % I put the hddot_max calculation only in those sections only because
            % otherwise when this function is used as an optimizer constraint, hddot_max
            % would be recomputed every step of the optimization.
            h(i) = kappa + (gamma + delta)*kappa_dot + hddot_max*(delta^2/2 + gamma*delta) ...
                + l_r*rho_r + l_v*rho_v;
        end
    else
        error('Unrecognized setting in CBF_obs');
    end
end
end