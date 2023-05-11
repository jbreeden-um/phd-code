function Plots
sim_lin15 = load('Runs/linear15.mat');
sim_lin30 = load('Runs/linear30.mat');
sim_lin45 = load('Runs/linear45.mat');
sim_lin60 = load('Runs/linear60.mat');
sim_nl15 = load('Runs/nonlin15.mat');
sim_nl30 = load('Runs/nonlin30.mat');
sim_nl45 = load('Runs/nonlin45.mat');
sim_nl60 = load('Runs/nonlin60.mat');
sim_nl90 = load('Runs/nonlin90.mat');
sim_nl120 = load('Runs/nonlin120.mat');
sim_nl180 = load('Runs/nonlin180.mat');
sim_nl240 = load('Runs/nonlin240.mat');
sim_nl300 = load('Runs/nonlin300.mat');
sim_nl360 = load('Runs/nonlin360.mat');
sim_nl420 = load('Runs/nonlin420.mat');
sim_nl480 = load('Runs/nonlin480.mat');
sim_nl540 = load('Runs/nonlin540.mat');
sim_p45 = load('Runs/planner45.mat');

figure(11); clf;
theta = linspace(0, 2*pi, 100);
cx = cos(theta);
cy = sin(theta);
rho = 1;
xc = get_center(0);
obs = obstacle_location(1,0)-xc(1:2); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r'); hold on;
obs = obstacle_location(2,0)-xc(1:2); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(3,0)-xc(1:2); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(4,0)-xc(1:2); fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
fill([-10 -10 10 10 -10], [6 0 0 6 6], 'r');
xlabel 'x_1 (km)'; ylabel 'x_2 (km)';

color1 = [0; 0.7; 0];
color2 = [0.2; 0.2; 1];
color3 = [0.9; 0; 0.9];
color4 = [0.4; 0.4; 0.4];
% p1 = plot(sim_lin15.x_lin(1,:)/1e3, sim_lin15.x_lin(2,:)/1e3, 'b', 'LineWidth', 2)
p2 = plot(sim_lin30.x_lin(1,:)/1e3, sim_lin30.x_lin(2,:)/1e3, '-', 'LineWidth', 2, 'Color', color1);
p3 = plot(sim_lin45.x_lin(1,:)/1e3, sim_lin45.x_lin(2,:)/1e3, '--', 'LineWidth', 2, 'Color', color1);
p5 = plot(sim_lin60.x_lin(1,:)/1e3, sim_lin60.x_lin(2,:)/1e3, ':', 'LineWidth', 3, 'Color', color1);
% p5 = plot(sim_nl15.x_lin(1,:)/1e3, sim_nl15.x_lin(2,:)/1e3, 'g', 'LineWidth', 2)
p6 = plot(sim_nl30.x_lin(1,:)/1e3, sim_nl30.x_lin(2,:)/1e3, '-', 'LineWidth', 2, 'Color', color2);
p7 = plot(sim_nl45.x_lin(1,:)/1e3, sim_nl45.x_lin(2,:)/1e3, '--', 'LineWidth', 2, 'Color', color2);
p8 = plot(sim_nl60.x_lin(1,:)/1e3, sim_nl60.x_lin(2,:)/1e3, ':', 'LineWidth', 3, 'Color', color2);
% p9 = plot(sim_nl90.x_lin(1,:)/1e3, sim_nl90.x_lin(2,:)/1e3, 'g', 'LineWidth', 2)
% p10 = plot(sim_nl120.x_lin(1,:)/1e3, sim_nl120.x_lin(2,:)/1e3, 'r', 'LineWidth', 2)
p11 = plot(sim_nl180.x_lin(1,:)/1e3, sim_nl180.x_lin(2,:)/1e3, '-', 'LineWidth', 2, 'Color', color3);
% p12 = plot(sim_nl240.x_lin(1,:)/1e3, sim_nl240.x_lin(2,:)/1e3, 'c', 'LineWidth', 2)
p13 = plot(sim_nl300.x_lin(1,:)/1e3, sim_nl300.x_lin(2,:)/1e3, '--', 'LineWidth', 2, 'Color', color3);
% p14 = plot(sim_nl360.x_lin(1,:)/1e3, sim_nl360.x_lin(2,:)/1e3, 'g', 'LineWidth', 2)
p15 = plot(sim_nl420.x_lin(1,:)/1e3, sim_nl420.x_lin(2,:)/1e3, ':', 'LineWidth', 3, 'Color', color3);
% p16 = plot(sim_nl480.x_lin(1,:)/1e3, sim_nl480.x_lin(2,:)/1e3, 'g', 'LineWidth', 2)
% p17 = plot(sim_nl540.x_lin(1,:)/1e3, sim_nl540.x_lin(2,:)/1e3, 'g', 'LineWidth', 2)
p18 = plot(sim_p45.x_lin(1,:)/1e3, sim_p45.x_lin(2,:)/1e3, '--', 'LineWidth', 2, 'Color', color4);

plot(sim_lin15.x0(1)/1e3, sim_lin15.x0(2)/1e3, 'ko', 'MarkerFaceColor', 'k');
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');

ax = gca;
legend(ax,[p2 p3 p5 p18], {'$$\psi_h$$, 30 sec', '$$\psi_h$$, 45 sec', '$$\psi_h$$, 60 sec', 'Plan, 45 sec'}, ...
    'FontSize', 13, 'Location', 'SouthEast', 'interpreter', 'latex', 'Position', [0.575257687635491 0.131746031746032 0.253313740935938 0.211190472103301]);
ax2 = axes('position',get(gca,'position'),'visible','off');
l2 = legend(ax2,[p6 p7 p8 p11 p13 p15], {'$$\psi_h^*$$, 30 sec', '$$\psi_h^*$$, 45 sec', '$$\psi_h^*$$, 60 sec', ...
    '$$\psi_h^*$$, 180 sec', '$$\psi_h^*$$, 300 sec', '$$\psi_h^*$$, 420 sec'}, ...
    'FontSize', 13, 'Location', 'SouthWest', 'interpreter', 'latex', 'Color', 'none');

set(l2, 'Position', [0.205357142857143 0.129365079365079 0.240648319691257 0.312619041488284]);
axis(ax,'equal');
axis(ax,[-8 8 -14 1]);
axis(ax2,'equal');
axis(ax2,[-8 8 -14 1]);

set(gcf, 'Position', [1972 271 560 420]);
text(0.25, -10.04, '$$x_0$$', 'interpreter', 'latex', 'FontSize', 14)
text(0.1, 0.55, '$$x_{target}$$', 'interpreter', 'latex', 'FontSize', 14)

%%
figure(12); clf;
i_end = 2000;
plot(sim_p45.t(1:i_end)/60, sim_p45.u(1,1:i_end), 'x', 'MarkerSize', 4, 'LineWidth', 1, 'Color', color4); hold on;
p1 = plot(sim_lin45.t(1:i_end)/60, sim_lin45.u(1,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', color1, 'Color', color1); hold on;
p2 = plot(sim_nl45.t(1:i_end)/60, sim_nl45.u(1,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', color2, 'Color', color2);
plot(sim_nl180.t(1:i_end)/60, sim_nl180.u(1,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', color3, 'Color', color3);
set(gca, 'Xlim', [0 30], 'YLim', [-12 30]);
set(gcf, 'Position', [2600 900 560 160]);
l = legend([p1 p2], {'$$\psi_h$$, 45 sec', '$$\psi_h^*$$, 45 sec'}, ...
    'FontSize', 12, 'Location', 'NorthEast', 'interpreter', 'latex');
xlabel 'Time (minutes)';
ylabel 'u_1 (m/s)';
set(l, 'Position', [0.677405256599883 0.617291670242945 0.214856648162022 0.268124996423721]);
% set(l, 'Position', [0.664113585070648 0.526666674713293 0.240648319691257 0.421249991953373]);


figure(13); clf;
p4 = plot(sim_p45.t(1:i_end)/60, sim_p45.u(2,1:i_end), 'x', 'MarkerSize', 4, 'LineWidth', 1, 'Color', color4); hold on;
plot(sim_lin45.t(1:i_end)/60, sim_lin45.u(2,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [0;0.7;0], 'Color', color1); hold on;
plot(sim_nl45.t(1:i_end)/60, sim_nl45.u(2,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [0;0;1], 'Color', color2);
p3 = plot(sim_nl180.t(1:i_end)/60, sim_nl180.u(2,1:i_end), 's', 'MarkerSize', 4, 'MarkerFaceColor', [1;0;1], 'Color', color3);
set(gca, 'Xlim', [0 30], 'YLim', [-30 15]);
set(gcf, 'Position', [2600 600 560 160]);
l = legend([p3 p4], {'$$\psi_h^*$$, 180 sec', 'Plan, 45 sec'}, ...
    'FontSize', 12, 'Location', 'SouthEast', 'interpreter', 'latex');
xlabel 'Time (minutes)';
ylabel 'u_2 (m/s)';
set(l, 'Position', [0.649642626366452 0.308333333333333 0.240833564109739 0.268124996423721]);

%%
figure(14); clf;
sim = sim_lin45;
p1a = plot(sim.t/60, sim.V, 'Color', color1); hold on;
p1b = plot(sim.t(sim.jumps)/60, sim.V(sim.jumps), '.', 'Color', color1, 'MarkerSize', 10);
sim = sim_nl45;
p2a = plot(sim.t/60, sim.V, 'Color', color2); hold on;
p2b = plot(sim.t(sim.jumps)/60, sim.V(sim.jumps), '.', 'Color', color2, 'MarkerSize', 10);
sim = sim_nl180;
p3a = plot(sim.t/60, sim.V, 'Color', color3); hold on;
p3b = plot(sim.t(sim.jumps)/60, sim.V(sim.jumps), '.', 'Color',  color3, 'MarkerSize', 10);
sim = sim_p45;
p4a = plot(sim.t/60, sim.V, 'Color', color4); hold on;
p4b = plot(sim.t(sim.jumps)/60, sim.V(sim.jumps), '.', 'Color',  color4, 'MarkerSize', 10);
p5 = plot(-1, 1, 'k.', 'MarkerSize', 10);
set(gcf, 'Position', [2600 300 560 180]);
set(gca, 'Xlim', [0 30]);
xlabel 'Time (minutes)';
ylabel 'V';
l = legend([p1a p2a p3a p4a p5], ...
    {'$$\psi_h$$, 45 sec, Flows', '$$\psi_h^*$$, 45 sec, Flows', ...
    '$$\psi_h^*$$, 18 sec, Flows', 'Plan, 45 sec, Flows', 'Jumps'}, ...
    'FontSize', 12, 'Location', 'NorthEast', 'interpreter', 'latex');
set(l, 'Position', [0.585686504203057 0.392592600539878 0.330278144246991 0.566666658719380]);

%%
compute_lin = [sim_lin15.mean_compute;
             sim_lin30.mean_compute;
             sim_lin45.mean_compute;
             sim_lin60.mean_compute]
         
compute_nl = [sim_nl15.mean_compute;
            sim_nl30.mean_compute;
            sim_nl45.mean_compute;
            sim_nl60.mean_compute;
            sim_nl90.mean_compute;
            sim_nl120.mean_compute;
            sim_nl180.mean_compute;
            sim_nl240.mean_compute;
            sim_nl300.mean_compute;
            sim_nl360.mean_compute;
            sim_nl420.mean_compute;
            sim_nl480.mean_compute;
            sim_nl540.mean_compute]
        
compute_select = [sim_lin45.mean_compute; sim_nl45.mean_compute; sim_nl180.mean_compute]

%%
data = {sim_lin15, sim_lin30, sim_lin45, sim_lin60, sim_nl15, sim_nl30, sim_nl45, sim_nl60, ...
    sim_nl90, sim_nl120, sim_nl180, sim_nl240, sim_nl300, sim_nl360, sim_nl420, sim_nl480, sim_nl540};
t_crit = zeros(1, length(data));
modified_cost = zeros(1, length(data));
tol = 6.9e4;
for i=1:length(data)
    [~,index] = max(data{i}.V<=tol);
    if index==1
        t_crit(i) = Inf;
        costs = vecnorm(data{i}.u);
        modified_cost(i) = Inf; %sum(costs(~isnan(costs)));
    else
        t_crit(i) = data{i}.t(index);
        costs = vecnorm(data{i}.u(:,1:index));
        modified_cost(i) = sum(costs(~isnan(costs)));
    end
end
disp('t_crit = ');
fprintf('   %d\n',t_crit);
disp('modified_cost = ');
for i=1:length(data), k = floor(log10(modified_cost(i))); fprintf([blanks(5-k), ' %.2f\n'],modified_cost(i)); end
end