function MakePlots
global constants s_target
m1 = load('Results/Results_m1optim.mat');
m1l = load('Results/Results_m1local.mat');
m1g = load('Results/Results_m1global.mat');
m1go = load('Results/Results_m1gloopt.mat');
m0 = load('Results/Results_m0.mat');
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

margin4 = [data.margin_m4optim];
for i=(length(margin4)-1):-1:2
    if margin4(i)==0 && margin4(i-1)~=0 && margin4(i+1)==0, margin4(i)=1e-5; end
end
for i=2:(length(margin4)-1)
    if margin4(i)==0 && margin4(i+1)~=0 && margin4(i-1)==0, margin4(i)=1e-5; end
end

h1 = semilogy(t, [data.margin_m0],       '--', 'Color', color_nu0, 'LineWidth', 1); hold on;
h2 = semilogy(t, [data.margin_m1global], '--', 'Color', color_nu1, 'LineWidth', 1);
h3 = semilogy(t, [data.margin_m1local],   '-', 'Color', color_nu1, 'LineWidth', 1);
h4 = semilogy(t, [data.margin_m1gloopt],  '--', 'Color', color_nu2, 'LineWidth', 1);
h5 = semilogy(t, [data.margin_m1optim],   '-', 'Color', color_nu2, 'LineWidth', 1);
h6 = semilogy(t, [data.margin_m4global], '--', 'Color', color_eta, 'LineWidth', 1);
h7 = semilogy(t, margin4,                '-', 'Color', color_eta, 'LineWidth', 1);

xlabel('Time', 'FontSize', 12);
ylabel('Required Margins', 'FontSize', 12);
axis([0 50 1e-5 5e1])
l = legend([h1,h2,h3,h4,h5,h6,h7], ...
    {'$\nu_0^g(T)$', '$\nu_1^g(T)$', '$\nu_1^l(T,x)$', '$\nu_2^g(T)$', ...
    '$\nu_2^l(T,x)$', '$\nu_3^g(T)$', '$\nu_3^l(T,x)$'}, ...
    'Interpreter', 'latex', 'FontSize', 12);
set(l, 'Position', [0.76 0.225 0.17 0.71])

yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3])

%%
f = figure(9); clf;
set(f, 'Position', [2800 800 600 450]);
PlotCone(constants.s, constants.theta);
hold on;
[s1, s2, s3] = sphere(20);
surf(s1, s2, s3, 'FaceAlpha', 0);
% End the simulations that diverged early.
m0.H(m0.x(:,3)<-0.1) = NaN;      m0.x(m0.x(:,3)<0,:) = NaN;
m1g.H(m1g.x(:,3)<-0.1) = NaN;    m1g.x(m1g.x(:,3)<0,:) = NaN;
                                 m1go.x(m1go.x(:,3)<0,:) = NaN;
h1 = plot3(m0.x(:,1),   m0.x(:,2),   m0.x(:,3),   '--', 'Color', color_nu0, 'LineWidth', 3);
h2 = plot3(m1g.x(:,1),  m1g.x(:,2),  m1g.x(:,3),  '--', 'Color', color_nu1, 'LineWidth', 3);
h3 = plot3(m1l.x(:,1),  m1l.x(:,2),  m1l.x(:,3),   '-', 'Color', color_nu1, 'LineWidth', 3);
h4 = plot3(m1go.x(:,1), m1go.x(:,2), m1go.x(:,3), '--', 'Color', color_nu2, 'LineWidth', 3);
h5 = plot3(m1.x(:,1),   m1.x(:,2),   m1.x(:,3),    '-', 'Color', color_nu2, 'LineWidth', 3);
h6 = plot3(m4g.x(:,1),  m4g.x(:,2),  m4g.x(:,3),  '--', 'Color', color_eta, 'LineWidth', 3);
h7 = plot3(m4.x(:,1),   m4.x(:,2),   m4.x(:,3),    '-', 'Color', color_eta, 'LineWidth', 3);
plot3(s_target(1), s_target(2), s_target(3), 'go', 'MarkerFaceColor', [0;1;0]);
axis equal; view([0 0 1]);
l = legend([h1,h2,h3,h4,h5,h6,h7], ...
    {'Case $\phi_0^g$', 'Case $\phi_1^g$', 'Case $\phi_1^l$', 'Case $\phi_2^g$', ...
    'Case $\phi_2^l$', 'Case $\phi_3^g$', 'Case $\phi_3^l$'}, ...
    'Interpreter', 'latex', 'FontSize', 16);
set(l, 'Location', 'NorthWest');

%%
f = figure(10); clf;
set(f, 'Position', [1800 700 600 200]);

h1 = semilogy(t, m0.H, '--', 'Color', color_nu0, 'LineWidth', 4); hold on;
h2 = semilogy(t, m1g.H, '--', 'Color', color_nu1, 'LineWidth', 3);
h3 = semilogy(t, m1l.H, '-', 'Color', color_nu1, 'LineWidth', 1);
h4 = semilogy(t, m1go.H, '--', 'Color', color_nu2, 'LineWidth', 1.5);
h5 = semilogy(t, m1.H, '-', 'Color', color_nu2, 'LineWidth', 1);
h6 = semilogy(t, m4g.H, '--', 'Color', color_eta, 'LineWidth', 1);
h7 = semilogy(t, m4.H, '-', 'Color', color_eta, 'LineWidth', 1);

xlabel('Time', 'FontSize', 12);
ylabel('CBF Value', 'FontSize', 12);
axis([0 40 -2e0 -1e-4])
try
    l = legend([h1,h2,h3,h4,h5,h6,h7], ...
        {'Case $\phi_0^g\;$', 'Case $\phi_1^g$', 'Case $\phi_1^l$', 'Case $\phi_2^g$', ...
        'Case $\phi_2^l$', 'Case $\phi_3^g$', 'Case $\phi_3^l$'}, ...
        'Interpreter', 'latex', 'FontSize', 12, 'NumColumns', 2);
    set(l, 'Position', [0.54 0.49 0.3587 0.4129])
catch
    l = legend([h1,h2,h3,h4,h5,h6,h7], ...
        {'Case $\phi_0^g\;$', 'Case $\phi_1^g$', 'Case $\phi_1^l$', 'Case $\phi_2^g$', ...
        'Case $\phi_2^l$', 'Case $\phi_3^g$', 'Case $\phi_3^l$'}, ...
        'Interpreter', 'latex', 'FontSize', 10);
    set(l, 'Position', [0.75 0.40 0.16 0.6]);
end
% set(l, 'Position', [0.75 0.40 0.16 0.6]);
% set(l, 'Position', [0.54 0.49 0.3587 0.4129])
yticks([-1e3 -1e2 -1e1 -1e0 -1e-1 -1e-2 -1e-3 -1e-4])

%%
f = figure(11); clf;
set(f, 'Position', [1800 400 600 200]);

u_m0 = [m0.data.u];
u_m1g = [m1g.data.u];
u_m1l = [m1l.data.u];
u_m1go = [m1go.data.u];
u_m1 = [m1.data.u];
u_m4g = [m4g.data.u];
u_m4 = [m4.data.u];

% h1 = plot(t, u_m0(1,:), '--', 'Color', color_nu0, 'LineWidth', 1); hold on;
% h2 = plot(t, u_m1g(1,:), '--', 'Color', color_nu1, 'LineWidth', 1);
h3 = plot(t, u_m1l(1,:), '-', 'Color', color_nu1, 'LineWidth', 1); hold on;
% h4 = plot(t, u_m1go(1,:), '--', 'Color', color_nu2, 'LineWidth', 1);
h5 = plot(t, u_m1(1,:), '-', 'Color', color_nu2, 'LineWidth', 1);
h6 = plot(t, u_m4g(1,:), '--', 'Color', color_eta, 'LineWidth', 1);
h7 = plot(t, u_m4(1,:), '-', 'Color', color_eta, 'LineWidth', 1);
axis([0 40 -0.012 0.012]);
xlabel('Time', 'FontSize', 12);
ylabel('Control Input  u_1', 'FontSize', 12);
l = legend([h3,h5,h6,h7], ...
    {'Case $\phi_1^l$', 'Case $\phi_2^l$', 'Case $\phi_3^g$', 'Case $\phi_3^l$'}, ...
    'Interpreter', 'latex', 'FontSize', 10);
set(l, 'Position', [0.744 0.205 0.16 0.35]);

f = figure(12); clf;
set(f, 'Position', [1800 100 600 200]);
h3 = plot(t, u_m1l(2,:), '-', 'Color', color_nu1, 'LineWidth', 1); hold on;
h5 = plot(t, u_m1(2,:), '-', 'Color', color_nu2, 'LineWidth', 1);
h6 = plot(t, u_m4g(2,:), '--', 'Color', color_eta, 'LineWidth', 1);
h7 = plot(t, u_m4(2,:), '-', 'Color', color_eta, 'LineWidth', 1);
axis([0 40 -0.012 0.012]);
xlabel('Time', 'FontSize', 12);
ylabel('Control Input  u_2', 'FontSize', 12);
l = legend([h3,h5,h6,h7], ...
    {'Case $\phi_1^l$', 'Case $\phi_2^l$', 'Case $\phi_3^g$', 'Case $\phi_3^l$'}, ...
    'Interpreter', 'latex', 'FontSize', 10);
set(l, 'Position', [0.744 0.205 0.16 0.35]);

f = figure(13); clf;
set(f, 'Position', [2500 400 600 200]);
h3 = plot(t, u_m1l(3,:), '-', 'Color', color_nu1, 'LineWidth', 1); hold on;
h5 = plot(t, u_m1(3,:), '-', 'Color', color_nu2, 'LineWidth', 1);
h6 = plot(t, u_m4g(3,:), '--', 'Color', color_eta, 'LineWidth', 1);
h7 = plot(t, u_m4(3,:), '-', 'Color', color_eta, 'LineWidth', 1);
axis([0 40 -0.012 0.012]);
xlabel('Time', 'FontSize', 12);
ylabel('Control Input  u_3', 'FontSize', 12);
l = legend([h3,h5,h6,h7], ...
    {'Case $\phi_1^l$', 'Case $\phi_2^l$', 'Case $\phi_3^g$', 'Case $\phi_3^l$'}, ...
    'Interpreter', 'latex', 'FontSize', 10);
set(l, 'Position', [0.744 0.205 0.16 0.35]);

end

function arr = limit(arr, lim)
arr(arr < -lim) = -lim;
arr(arr > lim) = lim;
end

function print_plots
figure(8); print -depsc Space_MarginsM4;
figure(9); print -depsc Space_Trajectories;
figure(10); print -depsc Space_Barriers;
figure(11); print -depsc Space_Input1;
figure(12); print -depsc Space_Input2;
figure(13); print -depsc Space_Input3;
end