function [measurement, rho] = get_estimates(t, x)
persistent generator
if t==0
    generator = RandStream('mt19937ar','Seed',0);
end
V = compute_V(t,x);
rho_r = 1e-3*sqrt(V);
rho_v = 2e-5*sqrt(V);
dr = generator.randn(2,1);
dv = generator.randn(2,1);
if norm(dr) > 1, dr = dr / norm(dr); end
if norm(dv) > 1, dv = dv / norm(dv); end

measurement = x + [rho_r*dr; rho_v*dv];
rho = [rho_r; rho_v];
end