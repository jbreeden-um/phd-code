function Plots
sim1 = load('Runs/sim1.mat');
sim2 = load('Runs/sim2.mat');
sim3 = load('Runs/sim3.mat');
sim4 = load('Runs/sim4.mat');
sim5 = load('Runs/sim5.mat');
sim6 = load('Runs/sim6.mat');
sim7 = load('Runs/sim7.mat');
sim8 = load('Runs/sim8.mat');

figure(11); clf;
theta = linspace(0, 2*pi, 100);
cx = cos(theta);
cy = sin(theta);
rho = 1;
obs = obstacle_location(1,0); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r'); hold on;
obs = obstacle_location(2,0); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(3,0); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(4,0); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
fill([-10 -10 10 10 -10], [6 0 0 6 6], 'r');
xlabel 'x_1 (km)'; ylabel 'x_2 (km)';

% p1 = plot(sim1.x(1,:)/1e3, sim1.x(2,:)/1e3, 'b', 'LineWidth', 2)
p2 = plot(sim2.x(1,:)/1e3, sim2.x(2,:)/1e3, '-', 'LineWidth', 2, 'Color', [0;0;1])
p3 = plot(sim3.x(1,:)/1e3, sim3.x(2,:)/1e3, '--', 'LineWidth', 2, 'Color', [0;0;1])
% p4 = plot(sim4.x(1,:)/1e3, sim4.x(2,:)/1e3, 'b', 'LineWidth', 2)
% p5 = plot(sim5.x(1,:)/1e3, sim5.x(2,:)/1e3, 'b', 'LineWidth', 2)
p6 = plot(sim6.x(1,:)/1e3, sim6.x(2,:)/1e3, '-', 'LineWidth', 2, 'Color', [0;0.7;0])
p7 = plot(sim7.x(1,:)/1e3, sim7.x(2,:)/1e3, '--', 'LineWidth', 2, 'Color', [0;0.7;0])
% p8 = plot(sim8.x(1,:)/1e3, sim8.x(2,:)/1e3, 'b', 'LineWidth', 2)

plot(sim2.x0(1)/1e3, sim2.x0(2)/1e3, 'ko', 'MarkerFaceColor', 'k');
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');

legend([p2 p3 p6 p7], {'$$\psi_h$$, 15 sec', '$$\psi_h^*$$, 15 sec', '$$\psi_h$$, 30 sec', '$$\psi_h^*$$, 30 sec'}, ...
    'FontSize', 14, 'Location', 'SouthEast', 'interpreter', 'latex');
axis equal;
axis([-8 8 -14 1]);
text(0.25, -10.04, '$$x_0$$', 'interpreter', 'latex', 'FontSize', 14)
text(0.1, 0.55, '$$x_{target}$$', 'interpreter', 'latex', 'FontSize', 14)

figure(12); clf;
i_end = 2000;
% p2 = plot(sim2.t(1:i_end)/60, sim2.u(1,1:i_end), 'o', 'MarkerSize', 4, 'MarkerFaceColor', [0;0;1], 'Color', [0;0;1]); hold on;
p3 = plot(sim3.t(1:i_end)/60, sim3.u(1,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [0;0;1], 'Color', [0;0;1]); hold on;
% p6 = plot(sim6.t(1:i_end)/60, sim6.u(1,1:i_end), 'o', 'MarkerSize', 4, 'MarkerFaceColor', [0;0.7;0], 'Color', [0;0.7;0]);
p7 = plot(sim7.t(1:i_end)/60, sim7.u(1,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [0;0.7;0], 'Color', [0;0.7;0]);
set(gca, 'Xlim', [0 30], 'YLim', [-12 30]);
set(gcf, 'Position', [2600 900 560 160]);
legend([p3 p7], {'$$\psi_h^*$$, 15 sec', '$$\psi_h^*$$, 30 sec'}, ...
    'FontSize', 14, 'Location', 'NorthEast', 'interpreter', 'latex');
xlabel 'Time (minutes)';
ylabel 'u_1 (m/s)';

figure(13); clf;
plot(sim3.t(1:i_end)/60, sim3.u(2,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [0;0;1], 'Color', [0;0;1]); hold on;
plot(sim7.t(1:i_end)/60, sim7.u(2,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [0;0.7;0], 'Color', [0;0.7;0]);
set(gca, 'Xlim', [0 30], 'YLim', [-30 15]);
set(gcf, 'Position', [2600 600 560 160]);
xlabel 'Time (minutes)';
ylabel 'u_2 (m/s)';

figure(14); clf;
sim = sim2;
p2a = plot(sim.t/60, sim.V, 'Color', [0;0;0.5]); hold on;
p2b = plot(sim.t(sim.jumps)/60, sim.V(sim.jumps), '.', 'Color', [0;0;0.5], 'MarkerSize', 8);
sim = sim3;
p3a = plot(sim.t/60, sim.V, 'Color', [0.3;0.6;1]); hold on;
p3b = plot(sim.t(sim.jumps)/60, sim.V(sim.jumps), '.', 'Color', [0.3;0.6;1], 'MarkerSize', 8);
set(gcf, 'Position', [2600 300 560 180]);
set(gca, 'Xlim', [0 30]);
xlabel 'Time (minutes)';
ylabel 'V';
legend([p2a p2b p3a p3b], ...
    {'$$\psi_h$$, 15 sec, Flows', '$$\psi_h$$, 15 sec, Jumps', ...
    '$$\psi_h^*$$, 15 sec, Flows', '$$\psi_h^*$$, 15 sec, Jumps'}, ...
    'FontSize', 13, 'Location', 'NorthEast', 'interpreter', 'latex');

sim2.mean_compute
sim6.mean_compute

sim3.mean_compute
sim7.mean_compute

end