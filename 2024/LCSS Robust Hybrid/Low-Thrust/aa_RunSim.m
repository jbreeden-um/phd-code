% This is the main simulation file. Rune this file to produce the simulation results.

%%
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
clear F;

%% Initial Conditions
sim_dur = 24*3600*10;
[x0, oe] = get_center(0);
x0(5) = x0(5) + 1;
dT_measurement = 6*3600;
dt = 60;

global wg_rate wc_max
wg_rate = 0.02;
wc_max = 9.2e-6; % 1367 W/m2 * 2 m2 / 3e8 m/s

%% Simulate
N = sim_dur/dt + 120;
t = zeros(1, N)*NaN;
x = zeros(6, N)*NaN;
u = zeros(3, N)*NaN;
xhat = zeros(6, N)*NaN;
rhohat = zeros(2, N)*NaN;
hhat = zeros(20,N)*NaN;
compute = zeros(1,N)*NaN;
measurements = zeros(1, N, 'logical');
t(1) = 0;
x(:,1) = x0;
xhat(:,1) = x0;
rhohat(:,1) = [0;0];
curr_index = 1;
new_sim = 1;
t_curr = t(curr_index);
last_measurement = -Inf;
sim_start = tic;
while t_curr < sim_dur
    % Get estimates
    if t(curr_index) - last_measurement >= dT_measurement
        [xhat(:,curr_index+1), rhohat(:,curr_index+1)] = get_estimates(t(curr_index), x(:,curr_index));
        measurements(curr_index+1) = 1;
        t(curr_index+1) = t(curr_index);
        x(:,curr_index+1) = x(:,curr_index);
        last_measurement = t(curr_index);
        curr_index = curr_index + 1;
    end
    
    % Get control
    t_curr = t(curr_index);
    x_curr = x(:,curr_index);
    xhat_curr = xhat(:,curr_index);
    rhohat_curr = rhohat(:,curr_index);
    [u_curr, compute(curr_index), h_curr] = CalculateU(t_curr,xhat_curr,rhohat_curr);
    u(:,curr_index) = u_curr;
    hhat(:,curr_index) = h_curr;
    
    % Propagate
    t(curr_index+1) = t_curr + dt;
    x(:,curr_index+1) = UpdateX(t_curr,x_curr,u_curr,dt,1);
    [xhat(:,curr_index+1), rhohat(:,curr_index+1)] = UpdateEstimates(t_curr,xhat_curr,u_curr,rhohat_curr,dt);
            
    curr_index = curr_index + 1;
    waitbar(t_curr/sim_dur);
end
toc(sim_start);

%% Plots
if length(t) > N
    warning('More measurements were taken than expected. Consider pre-allocating more space');
    N = length(t);
end
total_cost = sum(vecnorm(u(:,~isnan(u(1,:)))))*dt;
mean_compute = mean(compute(~isnan(compute)));

[~,i_end] = max(isnan(t));
i_end = i_end-1;
if i_end < 1, i_end = length(t); end

x_lin = convert_to_hcw(t,x);
xhat_lin = convert_to_hcw(t,xhat);
figure(1); clf;
plot(t, x_lin(1,:), 'r'); hold on;
plot(t, xhat_lin(1,:), 'r--');
plot(t, x_lin(2,:), 'g');
plot(t, xhat_lin(2,:), 'g--');
plot(t, x_lin(3,:), 'b');
plot(t, xhat_lin(3,:), 'b--');
ylabel 'Position';

figure(2); clf
plot(t, x_lin(4,:), 'r'); hold on;
plot(t, xhat_lin(4,:), 'r--');
plot(t, x_lin(5,:), 'g');
plot(t, xhat_lin(5,:), 'g--');
plot(t, x_lin(6,:), 'b');
plot(t, xhat_lin(6,:), 'b--'); 
ylabel 'Velocity';

figure(3); clf;
u_lin = convert_to_hcw(t(1:i_end),u(:,1:i_end),1);
plot(t(1:i_end), u_lin(1,:), 'r'); hold on;
plot(t(1:i_end), u_lin(2,:), 'g');
plot(t(1:i_end), u_lin(3,:), 'b'); 
ylabel 'Control Input';

figure(4); clf;
plot(xhat_lin(1,:), xhat_lin(2,:), 'm', 'LineWidth', 2); hold on;
plot(x_lin(1,:), x_lin(2,:), 'b', 'LineWidth', 2); hold on;
axis equal;

if new_sim
h = zeros(20, length(t))*NaN;
for i=1:i_end
    [~,~,~,h(:,i)] = CBF_icosa(t(i),x(:,i),rhohat(:,i),u(:,i));
    waitbar(i/i_end);
end
new_sim = 0;
end
figure(5); clf;
plot(t, h); ylabel 'Constraint'; hold on;
plot(t(measurements), h(1,measurements), 'bo');
plot(t, hhat, '--');
plot(t(measurements), hhat(1,measurements), 'ro');

figure(6); clf;
subplot(2,1,1); plot(t, rhohat(1,:)); ylabel '\rho_r'
subplot(2,1,2); plot(t, rhohat(2,:)); ylabel '\rho_v'