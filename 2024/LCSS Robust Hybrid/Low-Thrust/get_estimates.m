function [measurement, rho] = get_estimates(t, x)
persistent generator
if t==0
    generator = RandStream('mt19937ar','Seed',0);
end
rho_r = 10;
rho_v = 0.04;
dr = generator.randn(3,1);
dv = generator.randn(3,1);
if norm(dr) > 1, dr = dr / norm(dr); end
if norm(dv) > 1, dv = dv / norm(dv); end

measurement = x + [rho_r*dr; rho_v*dv];
rho = [rho_r; rho_v];
end