function out = H_func(x)
global constants
theta = constants.theta;
s = constants.s;
mu = constants.mu;
r = x(1:3);
w = x(4:6);
out = -cos(theta) + dot(r,s) + absSq(dot(s, cross(w, r)))/(2*mu);
end