function [h,A,b,hhat] = CBF_dock(t,x,T,rhohat,u_norm,setting)
% Function for h_7
global psi_h_nh
persistent hddot_max

mu = 398600e9;
[xc, ~, oe] = get_center(t);
uc = two_body(t,xc);
r = x(1:2) - xc(1:2);
v = x(3:4) - xc(3:4);
v_ref = sqrt(mu/oe.a);

vec = xc(3:4)/v_ref;
kappa = dot(vec, r);
kappa_dot = dot(vec, v) + dot(uc(3:4)/v_ref, r);

gamma = 1;
h = kappa + gamma*kappa_dot;

% The following Lipschitz constants could be generalized to the trajectory more precisely,
% but instead we simply pad them by 2% to account for variation over the trajectory.
if nargout == 4 && nargin==6 && isequal(setting,'post')
    l_margin = 1;
    A = [];
    b = [];
else
    l_margin = 1.02;
end
l_r = l_margin*(norm(vec) + gamma*norm(uc(3:4))/v_ref);
l_v = l_margin*gamma*norm(vec);

hhat = h + l_r*rhohat(1) + l_v*rhohat(2);

if nargout == 2 || nargout == 3 || (nargin==6 && isequal(setting,'nonlin'))
    if nargin==5 || isequal(setting, 'lin')
        hddot_max = 0.000467*sqrt(compute_V(t,x));

        Lfh = vec';
        kappa_dot_unc = Lfh*v + dot(r,uc(3:4))/v_ref;
        rho_sim = UpdateEstimates_Flow(t,x,rhohat,T);

        b = -(kappa + hddot_max*(T^2/2 + gamma*T) + kappa_dot_unc*(T+gamma) + l_r*rho_sim(1) + l_v*rho_sim(2));
        A = Lfh*(T+gamma);
    elseif isequal(setting, 'nonlin')
        delta = T/psi_h_nh;

        [~, rhohat] = UpdateEstimates_Jump(t,x,rhohat,[u_norm;0]);
        [tspan, x_sim, rho_sim] = UpdateEstimates_Flow(t,x,rhohat,T,psi_h_nh+1);
        h = zeros(psi_h_nh,1);
        
        for i=1:psi_h_nh
            xc = get_center(tspan(i));
            uc = two_body(tspan(i),xc);
            vec = xc(3:4)/v_ref;
            rho_r = rho_sim(1,i+1);
            rho_v = rho_sim(2,i+1);
            kappa = dot(vec, x_sim(1:2,i) - xc(1:2));
            kappa_dot = dot(vec, x_sim(3:4,i) - xc(3:4)) + dot(uc(3:4)/v_ref, x_sim(1:2,i) - xc(1:2));
            
            % Note that for this formula to be accurate, hddot_max must be computed by
            % calling the function with the linear or robust setting every so often
            h(i) = kappa + (gamma + delta)*kappa_dot + hddot_max*(delta^2/2 + gamma*delta) ...
                + l_r*rho_r + l_v*rho_v;
        end
    else
        error('Unrecognized setting in CBF_dock');
    end

end
end