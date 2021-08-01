function [u, rvecg] = GroundMotion(x)
% Computes a path on the surface of the asteroid slightly ahead of the spacecraft, so that
% it keeps moving.
%   u is the control input of the ground vehicle, which acts like a single integrator
%   rvecg is the present location of the ground vehicle
persistent fr fr_p fr_t
if isempty(fr)
    % This file uses a spline model of Eros to generate a path around the asteriod.
    % The spline is not very good around the poles, but that is okay, since the simulation
    % does not go near the poles.
    Eros = load('InData/ErosFit.mat');
    fr = Eros.fr;
    fr_p = Eros.fr_p;
    fr_t = Eros.fr_t;
end
global ground_vec

rvec = x(1:3);
vvec = x(4:6);
s = x(7);

phi0 = -pi/12;
dphi = pi/6;
theta0 = pi/4;
dtheta = 6*pi;
rho = 2.5; % km separation from spacecraft
smax = 1;

phi = phi0 + dphi*s/smax;
theta = theta0 + dtheta*s/smax;
while theta > pi, theta = theta - 2*pi; end
rg = fnval(fr, [phi; theta]);
rvecg = rg*[cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
ground_vec = rvecg;

h = rho - norm(rvecg - rvec);
Lf = (rvecg - rvec)'*vvec/norm(rvecg - rvec);

dr_dp = rg*[-sin(phi)*cos(theta); -sin(phi)*sin(theta); cos(phi)] ...
    + fnval(fr_p, [phi; theta])*[cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
dr_dt = rg*[-cos(phi)*sin(theta); cos(phi)*cos(theta); 0] ...
    + fnval(fr_t, [phi; theta])*[cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
dp_ds = dphi/smax;
dt_ds = dtheta/smax;
dr_ds = dr_dp*dp_ds + dr_dt*dt_ds;

Lg = -(rvecg - rvec)'*dr_ds/norm(rvecg-rvec);

alpha = -h;

% Lf + Lg*u \leq alpha
u = (alpha - Lf)/Lg;
if Lg > 0
    disp('Nonconvexity detected in the ground target planning');
    u = abs(u);
elseif u < 0
    u = 0;
end
end