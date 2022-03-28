function make_plots
data1 = load('data_OPCBF_high_res_tm150');
data2 = load('data_ECBF');
data3 = load('data_no_action');

%% Control Input Plots
f = figure(2); clf;
plot(data1.t(1:end-1), vecnorm(data1.u)*1e3, 'b', 'LineWidth', 2.5); hold on;
plot(data2.t(1:end-1), vecnorm(data2.u)*1e3, 'g', 'LineWidth', 1.5);
xlabel 'Time (s)';
ylabel '||u|| (m/s^2)';
axis([0 2500 -4 55]);
set(f, 'Position', [2200 1100 560 180]);
legend({'OPCBF', 'ECBF'},'Location','NorthWest','FontSize',12);

%% Constraint Plots
f = figure(3); clf;
p1 = plot(data1.t_dist, data1.h_dist, 'b', 'LineWidth', 5); hold on;
p3 = plot(data3.t_dist, data3.h_dist, 'r', 'LineWidth', 2.5);
p2 = plot(data2.t_dist, data2.h_dist, 'g', 'LineWidth', 1.5);
xlabel 'Time (s)';
ylabel 'h (km)';
set(f, 'Position', [2200 800 560 180]);
axis([0 2500 -14e3 1e3]);
legend([p1 p2 p3],{'OPCBF', 'ECBF', 'None'},'Location','NorthWest','FontSize',12);