function out = H_func(x)
global constants
sigma = constants.sigma;
rho = constants.rho;
d = wrap_angle(x(3) - atan2(x(2), x(1)));
out = rho - sqrt(x(1)^2 + x(2)^2 - sigma*d^2);
end