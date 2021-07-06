data = load('Results/SimLanding.txt');
control = load('Results/ControlLanding.txt');

t = data(:,1);
r = data(:,2:4);
v = data(:,5:7);
% H = data(:,8);
H = control(:,2);
% u = data(:,9:11);
u = control(:,3:5);
sigma = data(:,12);
te = data(:,13);

rho = 476000;
altitude = vecnorm(r,2,2) - rho;
figure(1); clf;
plot(t, altitude/1e3)
xlabel 'Time (s)'; ylabel 'Altitude (km)'
axis tight
set(gcf, 'Position', [1850 950 560 200]);
set(gca, 'FontSize', 11)

figure(2); clf;
plot3(r(:,1)/1e3, r(:,2)/1e3, r(:,3)/1e3, 'LineWidth', 2);
[s1, s2, s3] = sphere(20);
hold on;
surf(s1*rho/1e3, s2*rho/1e3, s3*rho/1e3, 'FaceColor', [0.7;0.7;0.7], 'FaceAlpha', 1);
xlabel 'r_x (km)';
ylabel 'r_y (km)';
zlabel 'r_z (km)';
axis equal;
set(gcf, 'Position', [1280 900 560 420])
view(-37.5, 30)
set(gca, 'FontSize', 11)

figure(3); clf;
plot(t, u); hold on;
plot([0 t(end)],[-0.5 -0.5],'r--')
plot([0 t(end)],[0.5 0.5],'r--')
axis([0 t(end) -0.55 0.55]);
legend 'u_x' 'u_y' 'u_z'
l = legend; set(l, 'Orientation', 'horizontal');
xlabel 'Time (s)'; ylabel 'Control Input (m/s^2)';
set(gcf, 'Position', [1850 600 560 200]);
set(gca, 'FontSize', 11)

hdot_terminal = dot(r(end-1,:),v(end-1,:))/norm(r(end-1,:))