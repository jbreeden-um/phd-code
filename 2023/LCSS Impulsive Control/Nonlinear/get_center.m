function [x, theta] = get_center(t)
% Position of the target satellite
mu = 398600e9;
Re = 6378e3;
alt = 600e3;
n = sqrt(mu/(Re+alt)^3);
delta0 = 0*pi/4;
theta = n*t + delta0;
r = Re + alt;
v = sqrt(mu/r);
x = [r*cos(theta); r*sin(theta); -v*sin(theta); v*cos(theta)];
end