%% Non-Convex Demo Cases
try
    s = get_s(0);
    s(1);
catch
    ComputeConstants;
end
data_nom = load('Results/Nonconvex_Nominal.mat');
data_com = load('Results/Nonconvex_Comparison.mat');
data_two = load('Results/Nonconvex_Combined.mat');

color1 = [0; 0.6; 0];
color2 = [0.2; 0.2; 1];
color3 = [0.4; 0.4; 0];
color_w = [0.4; 0.4; 0.4];
color_u = [0.000, 0.447, 0.741; 0.850, 0.325, 0.098; 0.929, 0.694, 0.125; 0.494, 0.184, 0.556];    

t = data_nom.t;
wheel_limit = data_nom.constants.wheel_limit;
p1 = data_nom.constants.p1;
p2 = data_nom.constants.p2;
ctheta = data_nom.ctheta;
s_target = data_nom.s_target;

figure(1); clf;
x1 = RotQ(p1, data_nom.x(:,1:4)');
x2 = RotQ(p2, data_nom.x(:,1:4)');
x3 = RotQ(p1, data_com.x(:,1:4)');
x4 = RotQ(p2, data_com.x(:,1:4)');
x5 = RotQ(p1, data_two.x(:,1:4)');
x6 = RotQ(p2, data_two.x(:,1:4)');
az1 = atan2d(x1(2,:), x1(1,:));
az2 = atan2d(x2(2,:), x2(1,:));
az3 = atan2d(x3(2,:), x3(1,:));
az4 = atan2d(x4(2,:), x4(1,:));
az5 = atan2d(x5(2,:), x5(1,:));
az6 = atan2d(x6(2,:), x6(1,:));
az4(logical([0, abs(diff(az4)) > 1])) = NaN; % to prevent discontinuities
az6(logical([0, abs(diff(az6)) > 1])) = NaN;
el1 = asind(x1(3,:));
el2 = asind(x2(3,:));
el3 = asind(x3(3,:));
el4 = asind(x4(3,:));
el5 = asind(x5(3,:));
el6 = asind(x6(3,:));
[c1, c2, c3] = PlotCone(get_s(0), acos(ctheta));
indices = abs([c1(:), c2(:), c3(:)]*get_s(0) - ctheta) < 1e-5;
azc = atan2d(c2(indices), c1(indices));
elc = asind(c3(indices));
fill(azc,elc,'r','EdgeAlpha',0); hold on;
f1 = plot(az1, el1, 'LineWidth', 3, 'Color', color1);
f2 = plot(az2, el2, 'LineWidth', 3, 'Color', color2);
f3 = plot(az3, el3, '--', 'LineWidth', 2, 'Color', color1);
f4 = plot(az4, el4, '--', 'LineWidth', 2, 'Color', color2);
f5 = plot(az5, el5, ':', 'LineWidth', 2, 'Color', color1);
f6 = plot(az6, el6, ':', 'LineWidth', 2, 'Color', color2);
azs = atan2d(s_target(2), s_target(1));
els = asind(s_target(3));
plot(azs, els, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(az1(1), el1(1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
plot(az2(1), el2(1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
xlabel 'Azimuth in Inertial Frame (deg)'; ylabel 'Elevation in Inertial Frame (deg');
legend([f1,f2,f3,f4,f5,f6],...
    {'ZohCBF b_1','ZohCBF b_2','Comparison b_1','Comparison b_2','Combined b_1','Combined b_2'},...
    'FontSize',12);
axis equal;
axis([-180, 180, -90, 90]);
set(gcf, 'Position', [200 900 750 400]);
set(gca, 'FontSize', 12);

figure(2); clf;
f1 = plot(t, data_nom.H(1,:), 'Color', color1, 'LineWidth', 1); hold on;
f2 = plot(t, data_nom.H(2,:), 'Color', color2, 'LineWidth', 1);
f3 = plot(t, data_com.H(1,:), '--', 'Color', color1, 'LineWidth', 1);
f4 = plot(t, data_com.H(2,:), '--', 'Color', color2, 'LineWidth', 1);
f5 = plot(t, data_two.H(1,:), ':', 'Color', color1, 'LineWidth', 2);
f6 = plot(t, data_two.H(2,:), ':', 'Color', color2, 'LineWidth', 2);
xlabel 'Time (s)'; ylabel 'H';
% legend([f1,f2,f3,f4,f5,f6],{'CBF H_1','CBF H_2','L.F. H_1','L.F. H_2','CBF+L.F. H_1','CBF+L.F. H_2'},...
%     'Orientation','horizontal','Position',[0.145 0.29 0.75 0.15])
legend([f1,f2],{'H_1','H_2'},'Orientation','horizontal','Location','SouthEast')
axis([0 1100 -1.2 0]);
set(gcf, 'Position', [390 600 560 200]);
set(gca, 'FontSize', 12);

figure(3); clf;
f1 = plot(t, data_nom.h(3,:)*1e3, 'Color', color_w, 'LineWidth', 1); hold on;
f2 = plot(t, data_com.h(3,:)*1e3, '--', 'Color', color_w, 'LineWidth', 1);
f3 = plot(t, data_two.h(3,:)*1e3, ':', 'Color', color_w, 'LineWidth', 2);
xlabel 'Time (s)'; ylabel 'h_3 (kg-m^2/s^2)';
legend([f1,f2,f3],{'ZohCBF','Comparison','Combined'},'Position',[0.395 0.57 0.23 0.31]);
axis([0 1100 -0.06 0]);
set(gcf, 'Position', [390 300 560 200]);
set(gca, 'FontSize', 12);

for i=1:4
figure(3+i); clf;
plot(t, data_nom.u(i,:)*1e3, 'Color', color_u(i,:), 'LineWidth', 1); hold on;
plot(t, data_com.u(i,:)*1e3, '--', 'Color', color_u(i,:), 'LineWidth', 1);
plot(t, data_com.u(i,:)*1e3, ':', 'Color', color_u(i,:), 'LineWidth', 1);
plot([0 1100], [1 1]*wheel_limit*1e3, 'k--');
plot([0 1100], -[1 1]*wheel_limit*1e3, 'k--');
xlabel 'Time (s)'; ylabel(['u_' num2str(i) ' (mNm)']);
legend({'ZohCBF', 'Comparison', 'Combined'}, 'Position', [0.705 0.615 0.20 0.31]);
axis([0 350 -0.75 0.75]);
set(gcf, 'Position', [1000, 1100-300*(i-1), 650, 200]);
set(gca, 'FontSize', 12);
end

if 0
    figure(1); print -depsc AzElNonCon.eps;
    figure(2); print -depsc ConstraintQNonCon.eps;
    figure(3); print -depsc ConstraintENonCon.eps;
end