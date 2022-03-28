function [p, dtau, dt, dx] = path_func(tau,t,x)
p = propagate_orbit(x,tau-t);

if nargout > 1
    dtau = two_body(p,[0;0;0]);

    dt = -dtau;

    dx = zeros(6,6);
    delta = .01;
    for i=1:3
        x_new = x;
        x_new(i) = x_new(i) + delta;
        p_new = propagate_orbit(x_new,tau-t);
        dx(:,i) = (p_new-p)/delta;
    end
    delta = 0.001;
    for i=4:6
        x_new = x;
        x_new(i) = x_new(i) + delta;
        p_new = propagate_orbit(x_new,tau-t);
        dx(:,i) = (p_new-p)/delta;
    end
end

end

