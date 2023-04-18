% To end the simulation early, just close the waitbar. Error handling has been set up in
% such a way that it will plot the intermediary results.

%% Clean Up Workspace
clear CalculateU UpdateX ConstraintE ConstraintQ PlotOutputs
clear ComparisonControllerBarrier ComparisonControllerSMC ComparisonControllerMPC mpc_func3 mpc_solver mpc_PE
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% Gain Presets
% The following presents gains used in the ComparisonControllerBarrier
zoh_violation_demo = 0; % 0 or 1
nonconvex_demo = 1; % 0 or 1

%% Other Choices
global control_law use_cbf
control_law = 2; % 1 = PD, 2 = Logarithmic Barriers, 3 = SMC, 4 = MPC
use_cbf = 1; % whether or not we use the CBF and the QCQP
SimDur = 1100;

%% Set up Variables
% The following is used by the controller
global s_target k_compare
s_target = [0; -1/sqrt(2); -1/sqrt(2)];

constants = load('parameters.mat');
% The following constants are used by the functions h, hdot, phi
global O_Earth_from_Sol omega_Earth Z Jtot wheel_axis Jw ctheta
O_Earth_from_Sol = constants.O_Earth_from_Sol;
omega_Earth = constants.omega_Earth;
Z = constants.Z;
Jtot = constants.Jtot;
wheel_axis = constants.wheel_axis;
Jw = constants.Jw;
ctheta = constants.ctheta;

if zoh_violation_demo
    % A demo where the constraints get violated because the controller is
    % way too aggressive for the ZOH discretization used. This case is
    % obviously contrived and would be tuned away, but illustrates how a
    % check on the control input can be useful.
    use_cbf = 0;
    control_law = 2;
    k_compare = 1;
elseif nonconvex_demo
    % A demo where the shortest-path control law will fail because of the
    % nonconvexity of the problem.
    if use_cbf
        % a large but somewhat arbtirary gain
        k_compare = 0.04;
    else
        % largest gain such that the angular velocity constraint is not violated
        k_compare = 0.0104;
    end
else
    % Nominal demo case uses the largest k_compare such that the angular
    % velocity constraint is not violated
    k_compare = 0.0165; % Controller gain for the Lyapunov function based controller
end

%% Set Up Initial Conditions
q0 = [1/2; 1/2; 1/2; 1/2];
w0 = [0; 0; 0];
W0 = [0; 0; 0; 0];

if nonconvex_demo
    q0 = [1/2; 1/2; 0.9; 1/2];
    q0 = q0/norm(q0);
    SimDur = 1100;
    ctheta = cos(0.95);
end

x0 = [q0; w0; W0];

%% Run Simulation
opts = odeset('OutputFcn', @PlotOutputs);
opts.dt = constants.dt;
opts.RelTol = inf;
opts.debug = 1;
opts.compare = 1;
timer = tic;
[t, x, data] = ode01([0 SimDur], x0, opts);
toc(timer)

%% Plot Results
figure(2); clf;
x1 = RotQ(constants.p1, x(:,1:4)');
x2 = RotQ(constants.p2, x(:,1:4)');
plot(t, x1); hold on;
plot(t, x2, '--');
title 'Pointing Vector Trajectories';
xlabel 'Time (s)'; legend r_{1,x} r_{1,y} r_{1,z} r_{2,x} r_{2,y} r_{2,z} 

figure(7); clf;
az1 = atan2(x1(2,:), x1(1,:));
az2 = atan2(x2(2,:), x2(1,:));
el1 = asin(x1(3,:));
el2 = asin(x2(3,:));
[c1, c2, c3] = PlotCone(get_s(0), acos(ctheta));
indices = abs([c1(:), c2(:), c3(:)]*get_s(0) - ctheta) < 1e-5;
azc = atan2(c2(indices), c1(indices));
elc = asin(c3(indices));
g = fill(azc,elc,'r','EdgeAlpha',0); hold on;
plot(az1, el1, 'b', 'LineWidth', 2);
plot(az2, el2, 'c', 'LineWidth', 2);
azs = atan2(s_target(2), s_target(1));
els = asin(s_target(3));
plot(azs, els, 'go', 'MarkerFaceColor', 'g');
axis equal;
axis([-pi, pi, -pi/2, pi/2]);
converged = acos(x1'*s_target) <= deg2rad(0.1);
settling_time = min(t(converged))

figure(3); clf;
plot(t, x(:,5:7)); title 'Velocity Trajectory';
xlabel 'Time (s)'; ylabel 'v (rad/s)';

figure(4); clf;
h = [data.h];
H = [data.H];
plot(t, h([1,2],:), '--');
hold on; title 'Pointing Constraint Values';
plot(t, H([1,2],:));
xlabel 'Time (s)'; ylabel 'h';
legend h_1 h_2 H_1 H_2 Location SouthEast

figure(5); clf;
plot(t, h(3,:));
title 'System Energy Constraint';
xlabel 'Time (s)'; ylabel 'h (J)';

figure(6); clf;
u = [data.u];
plot(t, u);
hold on; title 'Control Input';
plot(t, ones(size(t))*constants.wheel_limit*[-1, 1], 'r--');
xlabel 'Time (s)'; ylabel 'u';
legend u_1 u_2 u_3 u_4 Location SouthEast

figure(7); clf;
compute = [data.compute];
plot(t(2:end),compute(2:end));
title 'Computation Time'; xlabel 'Time (s)'; ylabel 'Cost (s)';
[~,i_end] = max(converged); if i_end==1, i_end=length(t); end
i_end = i_end+500; if i_end>length(t), i_end=length(t); end
mean_compute = mean(compute(3:i_end))
max_compute = max(compute(3:i_end))

%% Clean up
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
clear F;