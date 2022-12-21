function [tau, h_of_tau, dx] = find_max(t,x)
% $$\boldsymbol{M}(t,x)$$
% Since this is a really simple example, assume there is only one maximum. If there were
% multiple maxima, this would probably return the first one anyways
global T
[tau, h_of_tau] = my_fmincon(@(z) -h_func(z,path_func(z,t,x)), t, t+T);
h_of_tau = -h_of_tau;

if nargout==3
    dx = zeros(1,4);
    delta = 0.01;
    for i=1:4
        x_new = x;
        x_new(i) = x_new(i) + delta;
        new_tau = my_fmincon(@(z) -h_func(z,path_func(z,t,x_new)), -0.2, tau+0.2);
        dx(i) = (new_tau - tau)/delta;
    end
end
end

function [xval, fval] = my_fmincon(func, t1, t2)
N = 10;
tol = 1e-5;
t = zeros(1,N);
f = zeros(1,N);

index = 1;
f(index) = func(t(index));
while t2 - t1 > tol
    t(:) = linspace(t1, t2, N);
    for i=1:N
        f(i) = func(t(i));
    end
    [~,index] = min(f);
    if index==1
        t1 = t(1);
        t2 = t(2);
    elseif index==N
        t1 = t(N-1);
        t2 = t(N);
    else
        t1 = t(index-1);
        t2 = t(index+1);
    end
end
xval = t(index);
fval = f(index);
end
