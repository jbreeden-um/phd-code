% To end the simulation early, just close the waitbar. Error handling has been set up in
% such a way that it will plot the intermediary results.

%% Clean Up Workspace
clear CalculateU UpdateX Constraint_Avoidance Constraint_AvoidanceNoLim Constraint_AvoidanceInfNorm Plot_Outputs
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% Set up Variables
SimDur = 100;
a_max = 0.04;
global sim_case
sim_case = 1; % Choose between {1, 2, 3}

%% Save Needed Data
save('InData/SimData.mat', 'a_max', 'SimDur');

%% Set Up Initial Conditions
r0 = [12; 8; 0];
v0 = [0; 0; 0];
s0 = 0;
x0 = [r0; v0; s0];

%% Run Simulation
opts = odeset('OutputFcn', @PlotOutputs);
func = @(t, x) UpdateX(t, x, CalculateU(t, x));
opts.dt = 1;
opts.RelTol = 1e-2;
opts.enforce_errors = 1;
opts.max_error_count = 20;
opts.debug = 1;
tic
[t, x, data] = ode11(func, [0 SimDur], x0, opts);
toc

%% Plot Results
figure(1); clf;
r = 10;
[s1, s2, s3] = sphere(20);
surf(s1*r, s2*r, s3*r, 'FaceColor', [.7; .7; .7], 'FaceAlpha', 0.6, 'EdgeAlpha', 0.7);
hold on;
plot3(x(:,1), x(:,2), x(:,3));
rc = [data.rc]';
axis equal;
xlabel x; ylabel y; zlabel z;

figure(2); clf;
drvec = x(:,1:3) - rc;
plot(t, vecnorm(drvec, 2, 2));
hold on; title 'Norm of r_c';
plot(t, 1*ones(size(t)));

figure(3); clf;
plot(t, x(:,1:3)); title 'State Trajectory';

f = figure(4); clf;
h = [data.h];
H = [data.H];
plot(t, h, '--');
hold on; title 'Constraint Values';
plot(t, H);
xlabel 'Time (t)'; ylabel 'h';
legend h_s 'H''' Location SouthEast
set(f, 'Position', [1200 200 560 250]);

a_max = 0.04;
f = figure(5); clf;
usol = [data.u];
plot(t, usol(1:3,:));
hold on; title 'Control Input';
plot(t, ones(size(t))*[-a_max, a_max], 'r--');
xlabel 'Time (t)'; ylabel 'u';
legend u_1 u_2 u_3 Location SouthEast
set(f, 'Position', [600 200 560 300]);

figure(6); clf;
plot(t, rc);
title 'Closest Point';

%% Clean up
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
clear F f