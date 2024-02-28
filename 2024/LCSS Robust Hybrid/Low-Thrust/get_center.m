function [x, oe_out] = get_center(t)
% Position of the target satellite
persistent oe M0
mu = 398600e9;
if t==0
    oe.a = 42164e3;
    oe.e = 0;
    oe.i = deg2rad(40);
    oe.Omega = 0;
    oe.omega = 0;
    oe.nu = 0;
    M0 = oe.nu;
end
n = sqrt(mu/oe.a^3);
M = M0 + n*t;
nu = M; % circular orbit
x = KeplerToCartesian(oe.a, oe.e, oe.i, oe.Omega, oe.omega, nu, mu);
oe_out = oe;
end