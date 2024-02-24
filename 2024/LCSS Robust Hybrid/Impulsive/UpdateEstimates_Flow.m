function [t, x, rho] = UpdateEstimates_Flow(t0, x0, rho0, T, N_to_return)
if nargout > 1
    [t,x] = UpdateX_Flow(t0,x0,T,N_to_return,0);
else
    t = [t0; t0+T];
    N_to_return = 2;
end
wc_max = 9.2e-6;
l_fr = 2.35e-6;
l_fv = 0;
A = [0, 1; l_fr, l_fv];
rho = zeros(2, N_to_return);
rho(:,1) = rho0;
delta_t = t(2)-t(1);
B1 = expm(A*delta_t);
B2 = inv(A)*(B1 - eye(2));
for i=2:N_to_return
    rho(:,i) = B1*rho(:,i-1) + B2*[0; wc_max];
end
if nargout == 1
    % In this case, just return the error, not the estimates.
    % This option is encoded to make nesting this inside an optimization loop easier.
    t = rho(:,2);
end
end