%% About Me
% Main simulation file. See the 2023 version of this file for some comments about
% feasibility of the converging part of the controller (not covered in the 2024 paper).
%
% Switch cases by adjusting TM
% 
% Note that the trajectory is plotted in a radius-other-normal frame, as computed in
% convert_to_vis

%%
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
clear F;

%% Timings
Ts = 3;
N_per_flow = 3;
TM = 60;
dT_compute = 6; % time delay
dT_after_control = TM - dT_compute;
N_before_control = dT_compute; % one sample per second
N_after_control = dT_after_control; % one sample per second
Ta = TM - 2*Ts;
Tm = 0;

%% Parameters explained in the prior 2023 paper
global psi_v_min_rate psi_v_rhs_max psi_v_tol_dt psi_v_tol_dT gamma_V psi_v_Nmax psi_h_nh
psi_v_tol_dt = 1e-6; % the amount of uncertainty we want the MPC predictions in psi_v to tolerate 
psi_v_tol_dT = 1e-4;
psi_v_Nmax = 12; % the maximum amount of computation time we will tolerate
psi_v_min_rate = -0.002; % =\gamma_1 in the paper
psi_v_rhs_max = 5e3;
gamma_V = 0.00002; % =(1-\gamma_2) in the paper
psi_h_nh = 10; % should be at least 2 or the code might break

%% Setup
mu = 398600e9;
[~, ~, oe] = get_center(0);
n = sqrt(mu/oe.a^3);
x0_hcw = [-7.4; -30e3; 0; -1];
x0 = x0_hcw + get_center(0) + [0; 0; [eye(2), zeros(2,1)]*cross([0;0;n], [x0_hcw(1:2); 0])];
sim_dur = 5000;

UpdateX_Jump(0,x0,[0;0],1); % Initialize jump disturbance model
for i=1:6, obstacle_location(i,0); end % Initialize obstacles in orbits
psi_v(Ts,0,x0,psi_v_tol_dt);

%% Run the simulation
N = sim_dur/Ts*N_per_flow + 120;
t = zeros(1, N)*NaN;
x = zeros(4, N)*NaN;
u = zeros(2, N)*NaN;
xhat = zeros(4, N)*NaN;
rhohat = zeros(2, N)*NaN;
compute = zeros(1,N)*NaN;
u_count = 1;
jumps = zeros(1, N, 'logical');
measurements = zeros(1, N, 'logical');
t(1) = 0;
x(:,1) = x0;
xhat(:,1) = x0;
rhohat(:,1) = [0;0];
t_curr = 0;
curr_index = 1;
total_cost = 0;
new_sim = 1;
sim_start = tic;
while t_curr < sim_dur
        % Get estimates
        t_curr = t(curr_index);
        x_curr = x(:,curr_index);
        [xhat_meas, rhohat_meas] = get_estimates(t_curr,x_curr);
        measurements(curr_index+1) = 1;

        % Propagate until control
        N_steps = N_before_control+1;
        [tt,xx] = UpdateX_Flow(t_curr, x_curr, dT_compute, N_steps, 1);
        t((curr_index+1):(curr_index+N_before_control+1)) = tt;
        x(:,(curr_index+1):(curr_index+N_before_control+1)) = xx;
        [~,yy,zz] = UpdateEstimates_Flow(t_curr, xhat_meas, rhohat_meas, dT_compute, N_steps);
        xhat(:,(curr_index+1):(curr_index+N_before_control+1)) = yy;
        rhohat(:,(curr_index+1):(curr_index+N_before_control+1)) = zz;
        curr_index = curr_index + N_steps;

        % Get new control
        t_curr = t(curr_index);
        x_control = x(:,curr_index);
        xhat_control = xhat(:,curr_index);
        rhohat_control = rhohat(:,curr_index);
        [u_curr, compute(u_count)] = CalculateU(t_curr, xhat_control, rhohat_control, TM, Ts, Tm, Ta);
        u_count = u_count+1;

        % Propagate until next measurement
        x_control = UpdateX_Jump(t_curr, x_control, u_curr, 1);
        [xhat_control, rhohat_control] = UpdateEstimates_Jump(t_curr, xhat_control, rhohat_control, u_curr);
        N_steps = N_after_control+1;
        [tt,xx] = UpdateX_Flow(t_curr, x_control, dT_after_control, N_steps, 1);
        jumps(curr_index + 1) = 1;
        t((curr_index+1):(curr_index+N_after_control+1)) = tt;
        x(:,(curr_index+1):(curr_index+N_after_control+1)) = xx;
        u(:,curr_index) = u_curr;
        [~,yy,zz] = UpdateEstimates_Flow(t_curr, xhat_control, rhohat_control, dT_after_control, N_steps);
        xhat(:,(curr_index+1):(curr_index+N_after_control+1)) = yy;
        rhohat(:,(curr_index+1):(curr_index+N_after_control+1)) = zz;
        curr_index = curr_index + N_steps;
        total_cost = total_cost + norm(u_curr);
    waitbar(t_curr/sim_dur);
end
toc(sim_start)

%% Plots
N = length(t);
total_cost;
mean_dwell = mean(diff(t(jumps)));
mean_compute = mean(compute(~isnan(compute)));

[~,i_end] = max(isnan(t));
i_end = i_end-1;
if i_end < 1, i_end = length(t); end

x_lin = convert_to_vis(t,x);
xhat_lin = convert_to_vis(t,xhat);
figure(1); clf;
plot(t, x_lin(1,:), 'b'); hold on;
plot(t, xhat_lin(1,:), 'b--');
plot(t, x_lin(2,:), 'r');
plot(t, xhat_lin(2,:), 'r--');
ylabel 'Position';

figure(2); clf
plot(t, x_lin(3,:), 'b'); hold on;
plot(t, xhat_lin(3,:), 'b--');
plot(t, x_lin(4,:), 'r');
plot(t, xhat_lin(4,:), 'r--'); 
ylabel 'Velocity';

V = zeros(1,N);
for i=1:N, V(i) = compute_V(t(i), x(:,i)); end
figure(3); clf;
plot(t, sqrt(V)); ylabel 'Sqrt Lyapunov Function';
hold on;
plot(t([jumps(2:end),jumps(1)]), sqrt(V([jumps(2:end),jumps(1)])), 'ro');
plot(t(jumps), sqrt(V(jumps)), 'bo');
legend 'Flow' 'Before Jump' 'After Jump';
[~,i] = max(V<=4.6);
if i==1, t_crit = Inf; V_end = V(i_end); else, t_crit = t(i); end

figure(4); clf;
plot(t(jumps), V([jumps(2:end),jumps(1)])./V(jumps), 'bx')
ylabel 'Lyapunov Contraction Ratio';

figure(5); clf;
theta = zeros(1, length(u));
for i=1:length(u)
    [~,theta(i)] = get_center(t(i));
end
u_lin = [u(1,:).*cos(theta) + u(2,:).*sin(theta); u(2,:).*cos(theta) - u(1,:).*sin(theta)];
plot(t(1:length(u)), u_lin, 'o'); ylabel 'Control Input';

figure(6); clf;
theta = linspace(0, 2*pi, 100);
cx = cos(theta);
cy = sin(theta);
rho = 1e3;
Pr = [eye(2), zeros(2)];
obs = convert_to_vis(0, obstacle_location(1,0)); fill(obs(1)+cx*rho, obs(2)+cy*rho, 'r'); hold on;
obs = convert_to_vis(0, obstacle_location(2,0)); fill(obs(1)+cx*rho, obs(2)+cy*rho, 'r');
obs = convert_to_vis(0, obstacle_location(3,0)); fill(obs(1)+cx*rho, obs(2)+cy*rho, 'r');
obs = convert_to_vis(0, obstacle_location(4,0)); fill(obs(1)+cx*rho, obs(2)+cy*rho, 'r');
obs = convert_to_vis(0, obstacle_location(5,0)); fill(obs(1)+cx*rho, obs(2)+cy*rho, 'r');
obs = convert_to_vis(0, obstacle_location(6,0)); fill(obs(1)+cx*rho, obs(2)+cy*rho, 'r');
xlabel 'x_1'; ylabel 'x_2';
plot(x_lin(1,:), x_lin(2,:), 'b', 'LineWidth', 2);
plot(xhat_lin(1,:), xhat_lin(2,:), 'm', 'LineWidth', 2);
axis equal;

if new_sim
h = zeros(7, length(t))*NaN;
hhat = zeros(7, length(t))*NaN;
for i=1:i_end
    h(1,i) = CBF_obs(t(i),x(:,i),TM,[0;0],[],1);
    h(2,i) = CBF_obs(t(i),x(:,i),TM,[0;0],[],2);
    h(3,i) = CBF_obs(t(i),x(:,i),TM,[0;0],[],3);
    h(4,i) = CBF_obs(t(i),x(:,i),TM,[0;0],[],4);
    h(5,i) = CBF_obs(t(i),x(:,i),TM,[0;0],[],5);
    h(6,i) = CBF_obs(t(i),x(:,i),TM,[0;0],[],6);
    h(7,i) = CBF_dock(t(i),x(:,i),TM,[0;0],[]);
    
    [~,~,~,hhat(1,i)] = CBF_obs(t(i),xhat(:,i),TM,rhohat(:,i),[],1,'post');
    [~,~,~,hhat(2,i)] = CBF_obs(t(i),xhat(:,i),TM,rhohat(:,i),[],2,'post');
    [~,~,~,hhat(3,i)] = CBF_obs(t(i),xhat(:,i),TM,rhohat(:,i),[],3,'post');
    [~,~,~,hhat(4,i)] = CBF_obs(t(i),xhat(:,i),TM,rhohat(:,i),[],4,'post');
    [~,~,~,hhat(5,i)] = CBF_obs(t(i),xhat(:,i),TM,rhohat(:,i),[],5,'post');
    [~,~,~,hhat(6,i)] = CBF_obs(t(i),xhat(:,i),TM,rhohat(:,i),[],6,'post');
    [~,~,~,hhat(7,i)] = CBF_dock(t(i),xhat(:,i),TM,rhohat(:,i),[],'post');
    waitbar(i/i_end);
end
new_sim = 0;
end
figure(7); clf;
plot(t, h(7,:), 'b'); ylabel 'Docking Constraint'; hold on;
plot(t(jumps), h(7,jumps), 'bo');
plot(t, hhat(7,:), 'b--');

figure(8); clf;
plot(t, h(1,:), 'r'); hold on;
plot(t, h(2,:), 'g');
plot(t, h(3,:), 'b');
plot(t, h(4,:), 'm');
plot(t, h(5,:), 'k');
plot(t, h(6,:), 'c'); ylabel 'Obstacle Constraints';
plot(t(jumps), h(1,jumps), 'ro');
plot(t(jumps), h(2,jumps), 'go');
plot(t(jumps), h(3,jumps), 'bo');
plot(t(jumps), h(4,jumps), 'mo');
plot(t(jumps), h(5,jumps), 'ko');
plot(t(jumps), h(6,jumps), 'co');
plot(t, hhat(1,:), 'r--');
plot(t, hhat(2,:), 'g--');
plot(t, hhat(3,:), 'b--');
plot(t, hhat(4,:), 'm--');
plot(t, hhat(5,:), 'k--');
plot(t, hhat(6,:), 'c--');

%%
figure(9); clf;
subplot(2,1,1); plot(t, rhohat(1,:)); ylabel '\rho_r'
subplot(2,1,2); plot(t, rhohat(2,:)); ylabel '\rho_v'