% Instructions: Run the file OtherNumbers.m. Then scroll to the bottom of OtherNumbers.m,
% and run the segment after the return statement. This populates the constants structure.
% Then select your method, and run this file.

% To end the simulation early, just close the waitbar. Error handling has been set up in
% such a way that it will plot the intermediary results.

% Note: the code will run faster if you comment out some of the unnecessary lines in
% Constraint.m

%% Clean Up Workspace
clear CalculateU UpdateX Constraint PlotOutputs
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% Choose a Method
global method
% method = 'm0';       % = \phi_0^g = Lemma 2
% method = 'm1local';  % = \phi_1^l = Theorem 1
% method = 'm1global'; % = \phi_1^g = Corollary 1
% method = 'm1optim';  % = \phi_2^l = Theorem 2
% method = 'm1gloopt'; % = \phi_2^g = Corollary 2
% method = 'm4optim';  % = \phi_3^l = Theorem 3
method = 'm4global'; % = \phi_3^g = Corollary 3

% Note some of the above methods will generate errors. That is perfectly normal. When
% errors occur, we choose to let the system either choose the best control input possible
% or let the system "coast" and see if it ever recovers.
% 
% The methods m1optim, m4optim, and m4global should generate no errors.
% For the other methods, the errors indicate feasibility issues, which are expected when
% the system is trying to achieve such a large margin. Guaranteed feasbility is an
% exciting area for future work.
% The error "Larger u may be required" means the QP was still successful despite the
% feasibility issue. Other warnings mean that the QP failed completely, so the control
% input was set to zero so the system would "coast".

%% Set up Variables
global s_target two_obstacles constants
constants.theta = pi/5;
constants.u_max = [0.01; 0.01; 0.01];
constants.w_max = [0.2; 0.2; 0.2];
s_target = [0; -1/sqrt(2); 1/sqrt(2)];
two_obstacles = 0;

SimDur = 40;

%% Set Up Initial Conditions
r0 = [1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];
w0 = [0; 0; 0];
x0 = [r0; w0];

%% Run Simulation
opts = odeset('OutputFcn', @PlotOutputs);
opts.dt = .1;
opts.debug = 1;
tic
[t, x, data] = ode01([0 SimDur], x0, opts);
toc

%% Plot Results
figure(2); clf;
plot(t, x(:,1:3)); title 'State Trajectory';

figure(3); clf;
plot(t, x(:,4:6)); title 'Velocity Trajectory';

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
plot(t, u);
hold on; title 'Control Input';
% plot(t, ones(size(t))*[-constants.u_max(1), constants.u_max(1)], 'r--');
xlabel 'Time (t)'; ylabel 'u';
legend u_1 u_2 Location SouthEast

figure(6); clf;
plot(t, [data.margin]);
xlabel 'Time (s)'; ylabel 'Assumed Margin (m/s)';

%% Clean up
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
% save(['Results/Results_' method]);