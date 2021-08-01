function [u, rvecg] = GroundMotion(x)
global ground_vec

rvec = x(1:3);
vvec = x(4:6);
s = x(7);

phi0 = -pi/12;
dphi = pi/6;
theta0 = pi/4;
dtheta = 6*pi;
rho = 2.5; % km separation from spacecraft
smax = 1; % Probably should normalize this to 1, but instead I'm scaling this to better suit numerics

phi = phi0 + dphi*s/smax;
theta = theta0 + dtheta*s/smax;
while theta > pi, theta = theta - 2*pi; end
rg = 10;
rvecg = rg*[cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
ground_vec = rvecg;

h = rho - norm(rvecg - rvec);
Lf = (rvecg - rvec)'*vvec/norm(rvecg - rvec);

dr_dp = rg*[-sin(phi)*cos(theta); -sin(phi)*sin(theta); cos(phi)];
dr_dt = rg*[-cos(phi)*sin(theta); cos(phi)*cos(theta); 0];
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