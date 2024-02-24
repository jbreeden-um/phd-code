function [x, theta, oe_out] = get_center(t)
% Position of the target satellite
persistent oe M0
mu = 398600e9;
if t==0
    Re = 6378e3;
    alt = 600e3;
    n = sqrt(mu/(Re+alt)^3);
    delta0 = 0;
    theta = n*t + delta0;
    r = Re + alt;
    v = sqrt(mu/r);
    v = v*1.05;
    x0 = [r*cos(theta); r*sin(theta); -v*sin(theta); v*cos(theta)];
    oe = CartesianToKepler([x0(1:2); 0], [x0(3:4); 0], mu);
    E0 = 2*atan(sqrt((1-oe.e)/(1+oe.e))*tan(oe.nu/2));
    M0 = E0 - oe.e*sin(E0);
end
n = sqrt(mu/oe.a^3);
M = M0 + n*t;
nu = convert_M_to_nu(M, oe.e);
y = KeplerToCartesian(oe.a, oe.e, oe.i, oe.Omega, oe.omega, nu, mu);
x = y([1;2;4;5]);
theta = atan2(y(2), y(1));
oe_out = oe;
end