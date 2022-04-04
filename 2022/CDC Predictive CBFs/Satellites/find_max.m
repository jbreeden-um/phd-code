function [tau, h_of_tau, dx] = find_max(t,x)
% $$\boldsymbol{M}(t,x)$$
% Since this is a really simple example, assume there is only one maximum. If there were
% multiple maxima, this would probably return the first one anyways
global T
[tau, h_of_tau] = fmincon(@(z) -h_func(z,path_func(z,t,x)), t, [], [], [], [], t, t+T, [], optimset('Display', 'off'));
h_of_tau = -h_of_tau;

if nargout==3
    dx = zeros(1,6);
    delta = 0.01;
    for i=1:3
        x_new = x;
        x_new(i) = x_new(i) + delta;
        new_tau = fmincon(@(z) -h_func(z,path_func(z,t,x_new)), t, [], [], [], [], t, t+T, [], optimset('Display', 'off'));
        dx(i) = (new_tau - tau)/delta;
    end
    delta = 0.001;
    for i=4:6
        x_new = x;
        x_new(i) = x_new(i) + delta;
        new_tau = fmincon(@(z) -h_func(z,path_func(z,t,x_new)), t, [], [], [], [], t, t+T, [], optimset('Display', 'off'));
        dx(i) = (new_tau - tau)/delta;
    end
end
end