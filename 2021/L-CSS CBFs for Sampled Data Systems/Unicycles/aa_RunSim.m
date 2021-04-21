% Instructions: Run the file OtherNumbers.m. Then scroll to the bottom of OtherNumbers.m,
% and run the segment after the return statement. This populates the constants structure.
% Then select your method, and run this file.

% To end the simulation early, just close the waitbar. Error handling has been set up in
% such a way that it will plot the intermediary results.

% Note: the code will run faster if you comment out some of the unnecessary lines in
% Constraint.m

%% Clean Up Workspace
clear CalculateU UpdateX Constraint Plot_Outputs
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% Choose a Method
global method two_obstacles
% method = 'm1local';  % = \phi_1^l = Theorem 1
% method = 'm1global'; % = \phi_1^g = Corollary 1
% method = 'm1optim';  % = \phi_2^l = Theorem 2
% method = 'm1gloopt'; % = \phi_2^g = Corollary 2
% method = 'm4optim';  % = \phi_3^l = Theorem 3
method = 'm4global'; % = \phi_3^g = Corollary 3

two_obstacles = 0; % 1 for yes, 0 for no

%% Set up Variables
global r_target separation constants
constants.rho = 10; % m
constants.u_max = [5; 0.25];
constants.u_min = [0; -0.25];
separation = 20.3;

SimDur = 100;

%% Set Up Initial Conditions
r0 = [-40; -30];
phi0 = pi/4;
x0 = [r0; phi0];
r_target = [20; 20];

%% Run Simulation
opts = odeset('OutputFcn', @PlotOutputs);
opts.dt = .1;
opts.RelTol = inf;
opts.debug = 1;
tic
[t, x, data] = ode01([0 SimDur], x0, opts);
toc

%% Plot Results
figure(1); clf;
xc = linspace(-constants.rho, constants.rho, 100);
yc = sqrt(constants.rho^2 - xc.^2);
fill([xc, fliplr(xc)], [yc, -yc], 'r', 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
hold on;
if two_obstacles
fill([xc, fliplr(xc)]-20.1/sqrt(2), [yc, -yc]+20.1/sqrt(2), 'r', 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
end
plot(x(:,1), x(:,2));
plot(r_target(1), r_target(2), 'go', 'MarkerFaceColor', [0;1;0]);
axis equal;
xlabel x; ylabel y;

figure(2); clf;
plot(t, x(:,3)); title 'Phi Trajectory';

figure(3); clf;
plot(t, x(:,1:2)); title 'State Trajectory';

try
figure(4); clf;
h = [data.h];
H = [data.H];
plot(t, h, '--');
hold on; title 'Constraint Values';
plot(t, H);
xlabel 'Time (t)'; ylabel 'h';
legend h 'H' Location SouthEast
catch
end

figure(5); clf;
u = [data.u];
plot(t, u(1:2,:));
hold on; title 'Control Input';
xlabel 'Time (t)'; ylabel 'u';
legend u_1 u_2 Location SouthEast

figure(6); clf;
plot(t, [data.margin]);
xlabel 'Time (s)'; ylabel 'Assumed Margin (m/s)';

%% Clean up
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F); clear F;
if two_obstacles==0
%     save(['Results/Results_' method]);
else
%     save(['Results/Two Obstacles/Results_' method]);
end