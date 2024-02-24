function [x, rho] = UpdateEstimates_Jump(t0, x0, rho0, u)
x = UpdateX_Jump(t0,x0,u,0);
rho = [rho0(1); rho0(2) + 0.05*norm(u)];
end