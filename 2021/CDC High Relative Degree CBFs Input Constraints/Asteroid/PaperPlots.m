%%
sim = load('OutData/Results.mat');
figure(1); clf;
Eros = load('InData/Eros_Shape.mat');
trisurf(Eros.plates+1, Eros.vertices(:,1), Eros.vertices(:,2), Eros.vertices(:,3), ...
    'FaceColor', [.7; .7; .7], 'FaceAlpha', 0.5, 'EdgeAlpha', 0.6);
hold on;
plot3(sim.x(:,1), sim.x(:,2), sim.x(:,3),'LineWidth',2);
axis equal;
xlabel 'x (km)'; ylabel 'y (km)'; zlabel 'z (km)';
view(179.3, 77.5);
set(gcf, 'Position', [200 350 1550 950]);
set(gca, 'FontSize', 20);

%%
figure(2); clf;
h = 1-vecnorm([sim.data.drc_alg]',2,2);
hdot = -dot([sim.data.drc_alg], sim.x(:,4:6)')'./(1-h);
H = h + max(hdot,0).^2/(sim.u_max*2);
plot(sim.t, h, '--'); hold on;
plot(sim.t, H);
title 'Constraint Values';
xlabel 'Time (s)'; ylabel 'h (km)'; legend max(\{h_a\}_i) max(\{H'\}_i) Location SouthEast
set(gcf, 'Position', [300 200 560 200]);
axis([0 600 -3 0]);

%%
figure(3); clf;
u = [sim.data.u];
plot(sim.t, u(1:3,:)*1e3);
hold on; title 'Control Input';
plot(sim.t, ones(size(sim.t))*[-sim.u_max, sim.u_max]*1e3, 'r--');
xlabel 'Time (s)'; ylabel 'u (m/s^2)';
set(gcf, 'Position', [1000 200 560 200]);
axis([0 600 -50 50]);
legend('u_x','u_y','u_z','Location',[0.754, 0.70, 0.115, 0.25]);

return
%% 
figure(1); print -depsc Draft15_Asteroid_PosPlot;
figure(2); print -depsc Draft15_Asteroid_ConstraintPlot;
figure(3); print -depsc Draft15_Asteroid_uPlot;