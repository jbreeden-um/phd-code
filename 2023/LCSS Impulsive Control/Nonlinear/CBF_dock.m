function [h,A,b] = CBF_dock(t,x,T,setting)
% Function for h_5
global psi_h_nh robust_margin_checking
persistent hddot_max

xc = get_center(t);
r = x(1:2) - xc(1:2);
v = x(3:4) - xc(3:4);

vec = xc(3:4)/norm(xc(3:4));
kappa = dot(vec, r);
kappa_dot = dot(vec, v);

gamma = 1;
h = kappa + gamma*kappa_dot;

if nargout > 1 || (nargin==4 && isequal(setting,'nonlin'))
    if nargin==3 || isequal(setting, 'lin')
        hddot_max = 0.0003444*sqrt(compute_V(t,x));

        Lfh = vec';
        kappa_dot_unc = Lfh*v;

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
            hddot_max = 0.0003444*sqrt(V);
        end
        
        for i=1:psi_h_nh
            xc = get_center(tspan(i));
            vec = xc(3:4)/norm(xc(3:4));
            kappa = dot(vec, x_sim(i,1:2)' - xc(1:2));
            kappa_dot = dot(vec, x_sim(i,3:4)' - xc(3:4));
            % Note that for this formula to be accurate, hddot_max must be computed by
            % calling the function with the linear or robust setting every so often
            h(i) = kappa + (gamma + delta)*kappa_dot + hddot_max*(delta^2/2 + gamma*delta);
        end
    else
        error('Unrecognized setting in CBF_dock');
    end

end
end