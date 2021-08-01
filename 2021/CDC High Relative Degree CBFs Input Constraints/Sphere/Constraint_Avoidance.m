function [A, b] = Constraint_Avoidance(x)
% Returns the A and b constraint matrices for ensuring avoidance is achieved.
% Results should be fed into QP as A \leq b
global is_error outdata
persistent a_max
if isempty(a_max)
    file = load('InData/SimData.mat');
    a_max = file.a_max;
end

rvec = x(1:3);
vvec = x(4:6);
r = norm(rvec);
v = norm(vvec);
rho = 11;
outdata.rc = rvec/norm(rvec)*10;

h = rho - r;
hdot = -dot(rvec, vvec)/norm(r);
LgLf = -rvec'/r;
Lf2 = (dot(rvec, vvec)^2 - r^2*v^2)/r^3;
outdata.h = h;

if h > 0
    is_error = 1;
end
if hdot < 0
    H = h;
    A = [];
    b = [];
    outdata.H = H;
    return
else
    H = h + hdot^2/(a_max*2);
    outdata.H = H;
end
if H > 0
    is_error = 1;
end

alpha = -H; % straight line of slope 1

A = hdot*LgLf/a_max;
b = alpha - hdot*(1+Lf2/a_max);
end