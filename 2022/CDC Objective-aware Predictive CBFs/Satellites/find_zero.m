function [eta, dx] = find_zero(tau,t,x)
eta = fzero(@(z) h_func(z,path_func(z,t,x)), (tau+t)/2, optimset('Display', 'off'));
if isnan(eta)
    % If we are having trouble finding eta, choose a new starting point
    eta = fzero(@(z) h_func(z,path_func(z,t,x)), tau-15, optimset('Display', 'off'));
end

if eta > tau || isnan(eta)
    % If we are having a lot of trouble finding eta, choose a smarter starting point
    z = linspace(t, tau, 100);
    hh = zeros(1,100);
    for i=1:100, hh(i) = h_func(z(i),path_func(z(i),t,x)); end
    [~, i] = min(abs(hh));
    start = z(i);
    eta = fzero(@(z) h_func(z,path_func(z,t,x)), start, optimset('Display', 'off'));
    if eta > tau || isnan(eta)
        error(['t = ' num2str(t) ': Something went wrong in R']);
    end
end
if eta < t
    eta = t;
    disp(['t = ' num2str(t) ': Something is wrong. This case should be accounted for in Hstar']);
end

if nargout==2
    [~, dh] = h_func(eta,path_func(eta,t,x));
    [~, dp_dtau, ~, dp_dx] = path_func(eta,t,x);
    dx = -dh*dp_dx/(dh*dp_dtau);
end
end