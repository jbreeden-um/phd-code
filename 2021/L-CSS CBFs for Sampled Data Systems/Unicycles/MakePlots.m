function MakePlots
m1 = load('Results/Results_m1optim.mat');
m1go = load('Results/Results_m1gloopt.mat');
m4 = load('Results/Results_m4optim.mat');
m4g = load('Results/Results_m4global.mat');

input = 'Results/Results_m4optim.mat';
load(input);

color_nu0 = [0;0;0];
color_nu1 = [0;0.6;0];
color_nu2 = [0;0;1];
color_eta = [0.9;0;0.5];

f = figure(8); clf;
set(f, 'Position', [1800 1000 600 200]);

h1 = semilogy(t, [data.margin_m1global], '--', 'Color', color_nu1, 'LineWidth', 1); hold on;
h2 = semilogy(t, [data.margin_m1local], '-', 'Color', color_nu1, 'LineWidth', 1);
h3 = semilogy(t, [data.margin_m1gloopt], '--', 'Color', color_nu2, 'LineWidth', 1);
h4 = semilogy(t, [data.margin_m1optim], '-', 'Color', color_nu2, 'LineWidth', 1);
h5 = semilogy(t, [data.margin_m4global], '--', 'Color', color_eta, 'LineWidth', 1);
h6 = semilogy(t, [data.margin_m4optim], '-', 'Color', color_eta, 'LineWidth', 1.5);

xlabel('Time', 'FontSize', 12);
ylabel('Required Margins', 'FontSize', 12);
axis([0 115 1e-3 1e3])
l = legend([h1,h2,h3,h4,h5,h6], ...
    {'$\nu_1^g(T)$', '$\nu_1^l(T,x)$', '$\nu_2^g(T)$', '$\nu_2^l(T,x)$', '$\nu_3^g(T)$', '$\nu_3^l(T,x)$'}, ...
    'Interpreter', 'latex', 'FontSize', 12);
set(l, 'Position', [0.81 0.295 0.17 0.63])

yticks([1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3])

%%
f = figure(9); clf;
set(f, 'Position', [2800 800 600 400]);
xc = linspace(-constants.rho, constants.rho, 100);
yc = sqrt(constants.rho^2 - xc.^2);
fill([xc, fliplr(xc)], [yc, -yc], 'r', 'FaceAlpha', 0.5, 'EdgeAlpha', 0); hold on;
h1 = plot(m1go.x(:,1), m1go.x(:,2), '--', 'Color', color_nu2, 'LineWidth', 1);
h2 = plot(m1.x(:,1), m1.x(:,2), '-', 'Color', color_nu2, 'LineWidth', 1);
h3 = plot(m4g.x(:,1), m4g.x(:,2), '--', 'Color', color_eta, 'LineWidth', 1);
h4 = plot(m4.x(:,1), m4.x(:,2), '-', 'Color', color_eta, 'LineWidth', 1);
plot(r_target(1), r_target(2), 'go', 'MarkerFaceColor', [0;1;0]);
axis equal;
xlabel('x_1', 'FontSize', 12);
ylabel('x_2', 'FontSize', 12);
l = legend([h1,h2,h3,h4], ...
    {'Case $\phi_2^g$', 'Case $\phi_2^l$', 'Case $\phi_3^g$', 'Case $\phi_3^l$'}, ...
    'Interpreter', 'latex', 'FontSize', 14);
set(l, 'Location', 'SouthEast');

%%
f = figure(10); clf;
set(f, 'Position', [1800 700 600 200]);

h1 = semilogy(t, m1go.H, '--', 'Color', color_nu2, 'LineWidth', 1); hold on;
h2 = semilogy(t, m1.H, '-', 'Color', color_nu2, 'LineWidth', 1);
h3 = semilogy(t, m4g.H, '--', 'Color', color_eta, 'LineWidth', 1);
h4 = semilogy(t, m4.H, '-', 'Color', color_eta, 'LineWidth', 1);

xlabel('Time', 'FontSize', 12);
ylabel('CBF Value', 'FontSize', 12);
axis([0 100 -1e2 -1e-4])
l = legend([h1,h2,h3,h4], ...
    {'Case $\phi_2^g$', 'Case $\phi_2^l$', 'Case $\phi_3^g$', 'Case $\phi_3^l$'}, ...
    'Interpreter', 'latex', 'FontSize', 12);
set(l, 'Location', 'NorthEast')
yticks([-1e3 -1e2 -1e1 -1e0 -1e-1 -1e-2 -1e-3 -1e-4])

%% 
m4_two = load('Results/Two Obstacles/Results_m4optim.mat');
m4g_two = load('Results/Two Obstacles/Results_m4global.mat');
m1_two = load('Results/Two Obstacles/Results_m1optim.mat');

f = figure(11); clf;
set(f, 'Position', [2800 300 600 400]);
xc = linspace(-constants.rho, constants.rho, 100);
yc = sqrt(constants.rho^2 - xc.^2);
fill([xc, fliplr(xc)], [yc, -yc], 'r', 'FaceAlpha', 0.5, 'EdgeAlpha', 0); hold on;
fill([xc, fliplr(xc)]-20.25/sqrt(2), [yc, -yc]+20.25/sqrt(2), 'r', 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
plot(r_target(1), r_target(2), 'go', 'MarkerFaceColor', [0;1;0]);axis equal;

h1 = plot(m1_two.x(:,1), m1_two.x(:,2), '-', 'Color', color_nu2, 'LineWidth', 1);
h2 = plot(m4g_two.x(:,1), m4g_two.x(:,2), '--', 'Color', color_eta, 'LineWidth', 1);
h3 = plot(m4_two.x(:,1), m4_two.x(:,2), '-', 'Color', color_eta, 'LineWidth', 1);
xlabel('x_1', 'FontSize', 12);
ylabel('x_2', 'FontSize', 12);
axis equal;
l = legend([h1,h2,h3], ...
    {'Case $\phi_2^l$', 'Case $\phi_3^g$', 'Case $\phi_3^l$'}, ...
    'Interpreter', 'latex', 'FontSize', 14);
set(l, 'Location', 'SouthEast');

%% 
f = figure(12); clf;
set(f, 'Position', [1800 400 600 200]);

u_m1g = [m1go.data.u];
u_m1 = [m1.data.u];
u_m4g = [m4g.data.u];
u_m4 = [m4.data.u];

h1 = plot(t, u_m1g(1,:), '--', 'Color', color_nu2, 'LineWidth', 1); hold on;
h2 = plot(t, u_m1(1,:), '-', 'Color', color_nu2, 'LineWidth', 1);
h3 = plot(t, u_m4g(1,:), '--', 'Color', color_eta, 'LineWidth', 1);
h4 = plot(t, u_m4(1,:), '-', 'Color', color_eta, 'LineWidth', 1);
axis([0 100 0 6]);
xlabel('Time', 'FontSize', 12);
ylabel('Control Input  u_1', 'FontSize', 12);
l = legend([h1,h2,h3,h4], ...
    {'Case $\phi_2^g$', 'Case $\phi_2^l$', 'Case $\phi_3^g$', 'Case $\phi_3^l$'}, ...
    'Interpreter', 'latex', 'FontSize', 12);
set(l, 'Location', 'NorthEast')

f = figure(13); clf;
set(f, 'Position', [1800 100 600 200]);
h1 = plot(t, u_m1g(2,:), '--', 'Color', color_nu2, 'LineWidth', 1); hold on;
h2 = plot(t, u_m1(2,:), '-', 'Color', color_nu2, 'LineWidth', 1);
h3 = plot(t, u_m4g(2,:), '--', 'Color', color_eta, 'LineWidth', 1);
h4 = plot(t, u_m4(2,:), '-', 'Color', color_eta, 'LineWidth', 1);
axis([0 100 -0.15 0.3]);
xlabel('Time', 'FontSize', 12);
ylabel('Control Input  u_2', 'FontSize', 12);
l = legend([h1,h2,h3,h4], ...
    {'Case $\phi_2^g$', 'Case $\phi_2^l$', 'Case $\phi_3^g$', 'Case $\phi_3^l$'}, ...
    'Interpreter', 'latex', 'FontSize', 12);
set(l, 'Location', 'NorthEast')

end

function print_plots
figure(8); print -depsc Uni_MarginsM4;
figure(9); print -depsc Uni_Trajectories;
figure(10); print -depsc Uni_Barriers;
figure(11); print -depsc Uni_TwoObstacles;
figure(12); print -depsc Uni_Input1;
figure(13); print -depsc Uni_Input2;
end
