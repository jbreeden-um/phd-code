% This is the file that generates ErosFit.mat
% This file generates a spline of Eros used to compute the path of the ground vehicle that
% the spacecraft is following.

% Load Data
load 'Eros_Shape.mat';

% Convert to Spherical Coordinates
r = vecnorm(vertices, 2, 2);
phi = asin(vertices(:,3) ./ r);
theta = atan2(vertices(:,2), vertices(:,1));

% Add a wrap around constraint to the theta values
tol = pi/2;
edge_lower = theta < -pi + tol;
edge_upper = theta > pi - tol;
r1 = [r; r(edge_lower); r(edge_upper)];
theta1 = [theta; theta(edge_lower)+2*pi; theta(edge_upper)-2*pi];
phi1 = [phi; phi(edge_lower); phi(edge_upper)];

% Add a wrap around constraint to the phi values
tol = pi/2;
edge_lower = phi < -pi/2 + tol;
edge_upper = phi > pi/2 - tol;
r2 = [r1; r(edge_lower); r(edge_upper)];
theta2 = [theta1; wrap_around(theta(edge_lower)+pi); wrap_around(theta(edge_upper)+pi)];
phi2 = [phi1; -pi-phi(edge_lower); pi-phi(edge_upper)];

% Construct Fit
f = fit([phi2, theta2], r2, 'cubicinterp');
[a1, a2] = meshgrid(-pi/2:.01:pi/2, -pi:.01:pi+.01);
a3 = f(a1, a2);
x = a3.*cos(a2).*cos(a1);
y = a3.*sin(a2).*cos(a1);
z = a3.*sin(a1);

% Make Plots
figure(1); clf;
surf(x, y, z, 'EdgeAlpha', 0.1);
axis equal;
hold on;
trisurf(plates+1, vertices(:,1), vertices(:,2), vertices(:,3), ...
    'FaceColor', [.7; .7; .7], 'FaceAlpha', 0.1)

%% Make a second fit
xv = -2.3:.01:2.3;
yv = -3.8:.01:3.8;
[xvals, yvals] = meshgrid(xv, yv);
zvals = f(xvals, yvals)';
pp = csapi({xv, yv}, zvals);

a1 = (-pi/2:.01:pi/2)'; a2 = (-pi:.01:pi+.01);
a3 = fnval(pp, {a1, a2});
x = a3.*cos(a2).*cos(a1);
y = a3.*sin(a2).*cos(a1);
z = a3.*sin(a1);

% Make Plots
figure(2); clf;
surf(x, y, z, 'EdgeAlpha', 0.1);
axis equal;
hold on;
trisurf(plates+1, vertices(:,1), vertices(:,2), vertices(:,3), ...
    'FaceColor', [.7; .7; .7], 'FaceAlpha', 0.1)

%% Make Outputs
fr = pp;
fr_p = fnder(pp, [1, 0]);
fr_t = fnder(pp, [0, 1]);
fr_pp = fnder(pp, [2, 0]);
fr_pt = fnder(pp, [1, 1]);
fr_tt = fnder(pp, [0, 2]);
fr_ppp = fnder(pp, [0, 3]);
fr_ppt = fnder(pp, [2, 1]);
fr_ptt = fnder(pp, [1, 2]);
fr_ttt = fnder(pp, [0, 3]);

%%
save('ErosFit.mat', 'fr', 'fr_p', 'fr_t', 'fr_pp', 'fr_pt', 'fr_tt', 'fr_ppp', 'fr_ppt', 'fr_ptt', 'fr_ttt');


function out = wrap_around(in)
out = in;
out(out < -pi) = out(out < -pi) + 2*pi;
out(out > pi) = out(out > pi) - 2*pi;
end