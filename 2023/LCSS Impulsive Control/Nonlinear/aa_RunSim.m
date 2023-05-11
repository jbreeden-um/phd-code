%%
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% Parameters we can set
global psi_v_min_rate psi_v_rhs_max psi_v_tol_dt psi_v_tol_dT gamma_V psi_v_Nmax
psi_v_tol_dt = 1e-6; % the amount of uncertainty we want the MPC predictions in psi_v to tolerate 
psi_v_tol_dT = 1e-4;
psi_v_Nmax = 12; % the maximum amount of computation time we will tolerate
psi_v_min_rate = -0.002; % =\gamma_1 in the paper
psi_v_rhs_max = 5e3;
gamma_V = 0.00002; % =(1-\gamma_2) in the paper

global use_psi_h_star psi_h_nh
use_psi_h_star = 1; % 0 or 1, whether to use the vector psi_h_star or the scalar psi_h 
psi_h_nh = 10;

global use_dock_constraint
use_dock_constraint = 1; % 0 or 1, determines whether h5 is enforced or not

global use_moving_obstacles
use_moving_obstacles = 0; % the simulation works well with time-varying obstacles, but it is harder to visualize

%% Timings
dt_flow = 3;
N_per_flow = 3;
dT_control = 15;
N_per_control = dT_control/dt_flow*N_per_flow;

%% Pre-Defined Cases
% Every plot shown in the paper can be derived solely by varying the following number
sim_case = 6; % only simulations 2, 3, 6, and 7 are included in the paper to avoid clutter

if sim_case > 0 && sim_case <= 8
    use_psi_h_star = bitget(sim_case, 1);
    use_dock_constraint = bitget(sim_case, 2);
    if bitget(sim_case, 3)
        dT_control = 30;
    else
        dT_control = 15;
    end
end

%% Setup
global P A_sys
mu = 398600e9;
Re = 6378e3;
alt = 600e3;
n = sqrt(mu/(Re+alt)^3);
A_sys = [0,     0, 1,    0;
         0,     0, 0,    1;
         3*n^2, 0, 0,    2*n;
         0,     0, -2*n, 0];
B_sys = [zeros(2); eye(2)];
R = eye(2)*100;
Q = eye(4);
[~, P, ~] = lqr(A_sys,B_sys,Q,R);

x0 = [-7.4; -10e3; 0; -1];
sim_dur = 5000;

%% Run the simulation
N = sim_dur/dt_flow*N_per_flow + 120;
t = zeros(1, N)*NaN;
x = zeros(4, N)*NaN;
u = zeros(2, N)*NaN;
compute = zeros(1,N)*NaN;
u_count = 1;
jumps = zeros(1, N, 'logical');
t(1) = 0;
x(:,1) = x0;
t_curr = 0;
curr_index = 1;
total_cost = 0;
new_sim = 1;
while t_curr < sim_dur
    x_control = x(:,curr_index);
%     try
        [u_curr, compute(u_count)] = CalculateU(t(curr_index), x_control, dT_control, dt_flow);
        u_count = u_count+1;
        if norm(u_curr) == 0
            N_steps = N_per_flow+1;
            [tt, xx] = UpdateX_Flow(t_curr, x_control, dt_flow, N_steps);
            t((curr_index+1):(curr_index+N_per_flow)) = tt(2:end);
            x(:,(curr_index+1):(curr_index+N_per_flow)) = xx(:,2:end);
            curr_index = curr_index + N_per_flow;
            t_curr = t(curr_index);
        else
            x_control = UpdateX_Jump(x_control, u_curr);
            N_steps = N_per_control+1;
            [tt,xx] = UpdateX_Flow(t_curr, x_control, dT_control, N_steps);
            jumps(curr_index + 1) = 1;
            t((curr_index+1):(curr_index+N_per_control+1)) = tt;
            x(:,(curr_index+1):(curr_index+N_per_control+1)) = xx;
            u(:,curr_index) = u_curr;
            curr_index = curr_index + N_steps;
            t_curr = t(curr_index);
            total_cost = total_cost + norm(u_curr);
        end
%     catch e
%         disp('Ending simulation early');
%         disp(e);
%         t_curr = sim_dur;
%     end
    waitbar(t_curr/sim_dur);
end

%% Plots
N = length(t);
total_cost
mean_dwell = mean(diff(t(jumps)))
mean_compute = mean(compute(~isnan(compute)))

[~,i_end] = max(isnan(t));
i_end = i_end-1;
if i_end < 1, i_end = length(t); end

figure(1); clf;
plot(t, x(1:2,:)); ylabel 'Position';

figure(2); clf
plot(t, x(3:4,:)); ylabel 'Velocity';

V = zeros(1,N);
for i=1:N, V(i) = x(:,i)'*P*x(:,i); end
figure(3); clf;
plot(t, V); ylabel 'Lyapunov Function';
hold on;
plot(t([jumps(2:end),jumps(1)]), V([jumps(2:end),jumps(1)]), 'ro');
plot(t(jumps), V(jumps), 'bo');
legend 'Flow' 'Before Jump' 'After Jump';
[~,i] = max(V<=P(1,1));
if i==1, t_crit = Inf; V_end = V(i_end), else, t_crit = t(i), end

figure(4); clf;
plot(t(jumps), V([jumps(2:end),jumps(1)])./V(jumps), 'bx')
ylabel 'Lyapunov Contraction Ratio';

figure(5); clf;
plot(t(1:length(u)), u, 'o'); ylabel 'Control Input';

figure(6); clf;
theta = linspace(0, 2*pi, 100);
cx = cos(theta);
cy = sin(theta);
rho = 1;
obs = obstacle_location(1,0); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r'); hold on;
obs = obstacle_location(2,0); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(3,0); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(4,0); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
xlabel 'x_1'; ylabel 'x_2';
plot(x(1,:)/1e3, x(2,:)/1e3, 'b', 'LineWidth', 2)

if new_sim
h = zeros(6, length(t))*NaN;
for i=1:i_end
    h(1,i) = CBF_obs(t(i),x(:,i),dT_control,1);
    h(2,i) = CBF_obs(t(i),x(:,i),dT_control,2);
    h(3,i) = CBF_obs(t(i),x(:,i),dT_control,3);
    h(4,i) = CBF_obs(t(i),x(:,i),dT_control,4);
    h(5,i) = CBF_dock(t(i),x(:,i),dT_control);
    waitbar(i/i_end);
end
new_sim = 0;
end
figure(7); clf;
plot(t, h(5,:), 'b'); ylabel 'Docking Constraint'; hold on;
plot(t(jumps), h(5,jumps), 'bo');

figure(8); clf;
plot(t, h(1,:), 'r'); hold on;
plot(t, h(2,:), 'g');
plot(t, h(3,:), 'b');
plot(t, h(4,:), 'm'); ylabel 'Obstacle Constraints';
plot(t(jumps), h(1,jumps), 'ro');
plot(t(jumps), h(2,jumps), 'go');
plot(t(jumps), h(3,jumps), 'bo');
plot(t(jumps), h(4,jumps), 'mo');