function [x, rho] = UpdateEstimates(t0, x0, u, rho0, dt)
global wg_rate wc_max
if nargout > 1
    x = UpdateX(t0,x0,u,dt,0);
end
u_norm = norm(u);
wg = wg_rate*u_norm;
w = wg + wc_max;

mu = 398600e9;
a = 42164e3;
l_fr = 2*mu/a^3; % really small number :)
l_fv = 0;
A = [0, 1; l_fr, l_fv];
B1 = expm(A*dt);
B2 = inv(A)*(B1 - eye(2));
rho = B1*rho0 + B2*[0; w];
if nargout == 1
    % In this case, just return the error, not the estimates
    x = rho;
end
end