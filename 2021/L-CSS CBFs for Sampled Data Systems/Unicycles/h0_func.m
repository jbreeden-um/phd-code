function out = h0_func(x)
global constants
rho = constants.rho;
r = x(1:2);
out = rho - norm(r);
end