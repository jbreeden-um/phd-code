function [eta, dx] = find_zero(tau,t,x)
eta = my_fzero(@(z) h_func(z,path_func(z,t,x)), t, tau);

if nargout==2
    [~, dh] = h_func(eta,path_func(eta,t,x));
    [~, dp_dtau, ~, dp_dx] = path_func(eta,t,x);
    dx = -dh*dp_dx/(dh*dp_dtau);
end
end

function out = my_fzero(func,t1,t2)
h1 = func(t1);
h2 = func(t2);
if h1 > 0
    out = t1;
    return
elseif h2 < 0
    error('Invalid bounds.');
end
tau = (t1+t2)/2;
h = func(tau);
count = 0;
while abs(h) > 1e-5
    count = count+1;
    if h < 0
        t1 = tau;
    else
        t2 = tau;
    end
    tau = (t1+t2)/2;
    h = func(tau);
    if count > 1000
        disp('Failure in my_fzero.');
        break;
    end
end
out = tau;
end