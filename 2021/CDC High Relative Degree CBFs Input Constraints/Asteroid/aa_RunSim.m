% To end the simulation early, just close the waitbar. Error handling has been set up in
% such a way that it will plot the intermediary results.

%% Clean Up Workspace
clear CalculateU UpdateX SavedData Constraint_Asteroid4c2 Plot_Outputs
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
initialize_data;

%% Set up Variables
SimDur = 600;
u_max = 0.04;

%% Save Needed Data
save('InData/SimData.mat', 'u_max', 'SimDur');

%% Set Up Initial Conditions
r0 = [10; 8; 0];
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
try
    Eros;
catch
    Eros = load('InData/Eros_Shape.mat');
end
trisurf(Eros.plates+1, Eros.vertices(:,1), Eros.vertices(:,2), Eros.vertices(:,3), ...
    'FaceColor', [.7; .7; .7], 'FaceAlpha', 0.5, 'EdgeAlpha', 0.6);
hold on;
plot3(x(:,1), x(:,2), x(:,3));
rc = [data.rc]';
axis equal;

figure(3); clf;
plot(t, x(:,1:3)); title 'State Trajectory';

f = figure(4); clf;
for i=1:length(t), if isempty(data(i).h_max), data(i).h_max = nan; end; end
for i=1:length(t), if isempty(data(i).H_max), data(i).H_max = nan; end; end
h = [data.h_max];
H = [data.H_max];
plot(t, h, '--');
hold on; title 'Constraint Values';
plot(t, H);
xlabel 'Time (t)'; ylabel 'h';
legend h_s 'H''' Location SouthEast
set(f, 'Position', [1200 200 560 250]);

figure(5); clf;
usol = [data.u];
plot(t, usol(1:3,:));
hold on; title 'Control Input';
plot(t, ones(size(t))*[-u_max, u_max], 'r--');

figure(7); clf;
plot(t, [data.n_avoid]);
title('Number of QP Constraints');

f = figure(8); clf;
h = 1-vecnorm([data.drc_alg]',2,2);
hdot = -dot([data.drc_alg], x(:,4:6)')'./(1-h);
H = h + max(hdot,0).^2/(u_max*2);
plot(t, h, '--'); hold on;
plot(t, H);
title 'Constraint Values';
xlabel 'Time (s)'; ylabel 'h'; legend max(h_a) max(H') Location SouthEast
set(f, 'Position', [1200 200 560 250]);

%% Clean up
clear f;
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);