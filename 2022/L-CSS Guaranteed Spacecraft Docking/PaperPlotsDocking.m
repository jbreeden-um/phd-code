data = load('Results/SimDocking.txt');
control = load('Results/ControlDocking.txt');

t = data(:,1);
r = data(:,2:3);
v = data(:,4:5);
H = control(:,2);
Hr = control(:,3);
Hl = control(:,4);
Hv = control(:,5);
u = control(:,6:7);

Delta = 0.03;
offset = -r(:,2);
figure(1); clf;
plot(t, offset)
xlabel 'Time (s)'; ylabel 'In-Track Distance (m)'
axis tight
set(gcf, 'Position', [1850 950 560 200]);
set(gca, 'FontSize', 11)

figure(2); clf;
plot(t, u); hold on;
plot([0 t(end)],[-0.082 -0.082],'r--')
plot([0 t(end)],[0.082 0.082],'r--')
axis([0 t(end) -0.1 0.1]);
legend 'u_1' 'u_2'
l = legend; set(l, 'Orientation', 'horizontal');
xlabel 'Time (s)'; ylabel 'Control Input (m/s^2)';
set(gcf, 'Position', [1850 600 560 200]);
set(gca, 'FontSize', 11)

figure(3); clf;
plot(t, r(:,1)); hold on;
plot(t, Delta*ones(size(t)), 'r--');
plot(t, -Delta*ones(size(t)), 'r--');
xlabel 'Time (s)'; ylabel 'Cross-Track  x_1 (m)'
axis tight
axis([0 t(end) -0.035 0.035]);
set(gcf, 'Position', [1850 250 560 200]);
set(gca, 'FontSize', 11)

figure(4); clf;
n = sqrt(398600e9/(6378e3 + 400e3)^3);
r0 = (6378e3 + 400e3)*[cos(n*t), sin(n*t)];
plot(r0(:,1), r0(:,2), 'b', 'LineWidth', 2); hold on;
plot(r0(:,1)+r(:,1), r0(:,2)+r(:,2), 'g', 'LineWidth', 2);
xlabel 'r_x (m)';
ylabel 'r_y (m)';
axis equal;

figure(5); clf;
plot(t, [H, Hr, Hl, Hv]);
xlabel 'Time (s)'; ylabel 'H';

figure(6); clf;
plot(t, r);
xlabel 'Time (s)';
ylabel 'Position (m)';

hdot_terminal = v(end-1,2)