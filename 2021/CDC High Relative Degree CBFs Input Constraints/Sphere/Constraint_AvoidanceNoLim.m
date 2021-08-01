function [A, b] = Constraint_AvoidanceNoLim(x)
% Returns the A and b constraint matrices for ensuring avoidance is achieved.
% Results should be fed into QP as A \leq b
global is_error outdata

rvec = x(1:3);
vvec = x(4:6);
r = norm(rvec);
v = norm(vvec);
rho = 11;

h = rho - r;
hdot = -dot(rvec, vvec)/norm(r);
LgLf = -rvec'/r;
Lf2 = (dot(rvec, vvec)^2 - r^2*v^2)/r^3;

f = @(t) atan(t) + pi;
df = @(t) 1/(1+t^2);

h_r = f(-hdot)*h;
Lfh_r = -df(-hdot)*h*Lf2 + f(-hdot)*hdot;
Lgh_r = -df(-hdot)*h*LgLf;

alpha = -h_r; % straight line of slope 1

A = Lgh_r;
b = alpha - Lfh_r;

outdata.h = h;
outdata.H = h_r;
outdata.rc = rvec/r*10;
if h_r > 0
    is_error = 1;
end
end