function out = h0dot_func(x)
global constants
s = constants.s;
r = x(1:3);
w = x(4:6);
out = dot(cross(w,r),s);
end