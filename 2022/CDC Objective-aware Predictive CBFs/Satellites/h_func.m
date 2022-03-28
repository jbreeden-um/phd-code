function [h, dt, dx] = h_func(t,x)
h = h_func_real(t,x);

if nargout > 1
    delta = 1;
    dt = (h_func_real(t+delta,x) - h)/delta;
    
    dx = zeros(1,6);
    delta = 0.01;
    for i=1:3
        x_new = x;
        x_new(i) = x_new(i) + delta;
        dx(i) = (h_func_real(t,x_new) - h)/delta;
    end
    delta = 0.001;
    for i=4:6
        x_new = x;
        x_new(i) = x_new(i) + delta;
        dx(i) = (h_func_real(t,x_new) - h)/delta;
    end
end
end

function out = h_func_real(t,x)
rho = 1;
out = rho - norm(x(1:3) - second_orbit(t));
end

function out = second_orbit(t)
y = propagate_debris(t);
out = y(1:3);
end

