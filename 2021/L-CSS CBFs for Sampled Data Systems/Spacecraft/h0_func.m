function out = h0_func(x)
global constants
theta = constants.theta;
s = constants.s;
r = x(1:3);
out = -cos(theta) + dot(r,s);
end