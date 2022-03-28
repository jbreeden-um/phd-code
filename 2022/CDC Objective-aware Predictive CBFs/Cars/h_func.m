function [out, dx] = h_func(t,x)
out = h_func_real(x);

if nargout==2
    dx = zeros(1,4);
    delta = 0.01;
    for i=1:4
        x_new1 = x;
        x_new2 = x;
        x_new1(i) = x_new1(i) + delta;
        x_new2(i) = x_new2(i) - delta;
        dx(i) = (h_func_real(x_new1) - h_func_real(x_new2))/(2*delta);
    end
end
end

function out = h_func_real(x)
rho = 2;
out = rho - norm(lane1(x(1)) - lane2(x(3)));
end