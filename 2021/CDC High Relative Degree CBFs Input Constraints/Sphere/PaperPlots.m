sim1 = load('OutData/Results_Case1'); 
sim2 = load('OutData/Results_Case2');
sim3 = load('OutData/Results_Case3');

figure(1); clf;
r = 10;
[s1, s2, s3] = sphere(20);
surf(s1*r, s2*r, s3*r, 'FaceColor', [.7; .7; .7], 'FaceAlpha', 0.6, 'EdgeAlpha', 0.7);
hold on;
p1 = plot3(sim1.x(:,1), sim1.x(:,2), sim1.x(:,3), 'LineWidth', 2);
p2 = plot3(sim2.x(:,1), sim2.x(:,2), sim2.x(:,3), 'LineWidth', 2);
p3 = plot3(sim3.x(:,1), sim3.x(:,2), sim3.x(:,3), 'LineWidth', 2);
axis equal;
xlabel 'x (m)'; ylabel 'y (m)'; zlabel 'z (m)';
legend([p1,p2,p3],{'H', 'H''', 'H_o'},'FontSize',16);
view(-177.1, 71.6);
set(gcf, 'Position', [2100 431 1108 854]);
set(gca, 'FontSize', 16)

%%
figure(2); clf;
t1 = sim1.t;
t2 = sim2.t;
t3 = sim3.t;
h1 = [sim1.data.h];
H1 = [sim1.data.H];
h2 = [sim2.data.h];
H2 = [sim2.data.H];
h3 = [sim3.data.h];
H3 = [sim3.data.H];
plot(t1, h1, '--', 'LineWidth', 1); hold on; 
% plot(t2, h2, '--');
% plot(t3, h3, '--');
plot(t1, H1, 'LineWidth', 1);
% plot(t2, H2);
% plot(t3, H3);
title 'Constraint Values';
xlabel 'Time (s)'; ylabel 'h (m)';
legend h_a 'H' Location SouthEast
set(gcf, 'Position', [300 200 560 200]);
axis([0 100 -4 0]);

%%
a_max = 0.04;
u1 = [sim1.data.u];
u2 = [sim2.data.u];
u3 = [sim3.data.u];

figure(3); clf;
plot(t1, u1(1:3,:), 'LineWidth', 1);
hold on; title 'Control Input';
plot(t1, ones(size(t1))*[-a_max, a_max], 'r--', 'LineWidth', 1);
xlabel 'Time (s)'; ylabel 'u (m/s^2)';
set(gcf, 'Position', [1000 800 560 200]);
axis([0 100 -0.05 0.05])
legend('u_x','u_y','u_z','Location',[0.79, 0.63, 0.115, 0.35]);


figure(4); clf;
plot(t2, u2(1:3,:), 'LineWidth', 1);
hold on; title 'Control Input';
plot(t2, ones(size(t2))*[-a_max, a_max], 'r--', 'LineWidth', 1);
xlabel 'Time (s)'; ylabel 'u (m/s^2)';
set(gcf, 'Position', [1000 500 560 200]);
axis([0 100 -0.05 0.05])
legend('u_x','u_y','u_z','Location',[0.79, 0.63, 0.115, 0.35]);


figure(5); clf;
plot(t3, u3(1:3,:), 'LineWidth', 1);
hold on; title 'Control Input';
plot(t3, ones(size(t3))*[-a_max, a_max], 'r--', 'LineWidth', 1);
xlabel 'Time (s)'; ylabel 'u (m/s^2)';
set(gcf, 'Position', [1000 200 560 200]);
axis([0 100 -0.05 0.19])
legend('u_x','u_y','u_z','Location',[0.79, 0.63, 0.115, 0.35]);

return;
%%
figure(1); print -depsc Draft15_PosPlot
figure(2); print -depsc Draft15_ConstraintPlot
figure(3); print -depsc Draft15_uPlotInf
figure(4); print -depsc Draft15_uPlotNom
figure(5); print -depsc Draft15_uPlotCompare
