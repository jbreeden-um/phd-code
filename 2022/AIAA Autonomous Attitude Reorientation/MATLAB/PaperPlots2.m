%% Non-Convex Demo Cases
data_nom = load('Results/Nonconvex_CBF.mat');
data_com = load('Results/Nonconvex_Barrier.mat');
data_two = load('Results/Nonconvex_Combined.mat');

try
    cross(get_s(0),[1;0;0]);
catch
    ComputeConstants;
end

color1 = [0.2; 0.2; 1];
color2 = [0; 0.6; 0];
% color3 = [0.7; 0; 0];
% color4 = [0.5; 0.5; 0.5];
color5 = [1; 0.5; 0];
% color_w = [0.4; 0.4; 0.4];
color_u = [0.000, 0.447, 0.741; 0.850, 0.325, 0.098; 0.929, 0.694, 0.125; 0.494, 0.184, 0.556];    

t = data_nom.t;
wheel_limit = data_nom.constants.wheel_limit;
p1 = data_nom.constants.p1;
p2 = data_nom.constants.p2;
ctheta = data_nom.ctheta;
s_target = data_nom.s_target;

%%
figure(1); clf;
x1_1 = RotQ(p1, data_nom.x(:,1:4)');
x1_2 = RotQ(p2, data_nom.x(:,1:4)');
x2_1 = RotQ(p1, data_com.x(:,1:4)');
x2_2 = RotQ(p2, data_com.x(:,1:4)');
x5_1 = RotQ(p1, data_two.x(:,1:4)');
x5_2 = RotQ(p2, data_two.x(:,1:4)');
az1_1 = atan2d(x1_1(2,:), x1_1(1,:));
az1_2 = atan2d(x1_2(2,:), x1_2(1,:));
az2_1 = atan2d(x2_1(2,:), x2_1(1,:));
az2_2 = atan2d(x2_2(2,:), x2_2(1,:));
az5_1 = atan2d(x5_1(2,:), x5_1(1,:));
az5_2 = atan2d(x5_2(2,:), x5_2(1,:));
az2_2(logical([0, abs(diff(az2_2)) > 1])) = NaN; % to prevent discontinuities
az5_2(logical([0, abs(diff(az5_2)) > 1])) = NaN;
el1_1 = asind(x1_1(3,:));
el1_2 = asind(x1_2(3,:));
el2_1 = asind(x2_1(3,:));
el2_2 = asind(x2_2(3,:));
el5_1 = asind(x5_1(3,:));
el5_2 = asind(x5_2(3,:));
[c1, c2, c3] = PlotCone(get_s(0), acos(ctheta));
indices = abs([c1(:), c2(:), c3(:)]*get_s(0) - ctheta) < 1e-5;
azc = atan2d(c2(indices), c1(indices));
elc = asind(c3(indices));
fill(azc,elc,'r','EdgeAlpha',0); hold on;
f1_1 = plot(az1_1, el1_1, '-', 'LineWidth', 3, 'Color', color1);
f1_2 = plot(az1_2, el1_2, '--', 'LineWidth', 3, 'Color', color1);
f2_1 = plot(az2_1, el2_1, '-', 'LineWidth', 3, 'Color', color2);
f2_2 = plot(az2_2, el2_2, '--', 'LineWidth', 3, 'Color', color2);
f5_1 = plot(az5_1, el5_1, '-', 'LineWidth', 3, 'Color', color5);
f5_2 = plot(az5_2, el5_2, '--', 'LineWidth', 3, 'Color', color5);
azs = atan2d(s_target(2), s_target(1));
els = asind(s_target(3));
f7 = plot(azs, els, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
f8 = plot(az1_1(1), el1_1(1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
f9 = plot(az1_2(1), el1_2(1), 'ko', 'LineWidth', 3, 'MarkerSize', 8);
xlabel 'Azimuth in Inertial Frame (deg)'; ylabel 'Elevation in Inertial Frame (deg)';
legend([f1_1,f2_1,f5_1,f8,f9,f7],...
    {'ZohCBF','Log-B','Combined','Initial b_1','Initial b_2','Target for b_1'},...
    'FontSize',12,'Location','NorthEast');
axis equal;
axis([-180, 180, -90, 90]);
set(gcf, 'Position', [200 900 750 400]);
set(gca, 'FontSize', 14);

%%
figure(2); clf;
f1_1 = plot(t, data_nom.h(1,:), '-', 'Color', color1, 'LineWidth', 1); hold on;
f1_2 = plot(t, data_nom.h(2,:), '--', 'Color', color1, 'LineWidth', 1);
f2_1 = plot(t, data_com.h(1,:), '-', 'Color', color2, 'LineWidth', 1);
f2_2 = plot(t, data_com.h(2,:), '--', 'Color', color2, 'LineWidth', 1);
f5_1 = plot(t, data_two.h(1,:), '-', 'Color', color5, 'LineWidth', 1);
f5_2 = plot(t, data_two.h(2,:), '--', 'Color', color5, 'LineWidth', 1);
f1 = plot(0,0,'k-', 'LineWidth', 1);
f2 = plot(0,0,'k--', 'LineWidth', 1);
xlabel 'Time (s)'; ylabel '\kappa_{ }';
% legend([f1,f2,f3,f4,f5,f6],{'CBF H_1','CBF H_2','L.F. H_1','L.F. H_2','CBF+L.F. H_1','CBF+L.F. H_2'},...
%     'Orientation','horizontal','Position',[0.145 0.29 0.75 0.15])
legend([f1,f2],{'\kappa_1','\kappa_2'},'Orientation','horizontal','Location','SouthEast')
axis([0 1100 -1.2 0]);
set(gcf, 'Position', [390 600 560 180]);
ax2 = gca;
set(ax2, 'FontSize', 12);

figure(3); clf;
f1_1 = plot(t, data_nom.h(3,:)*1e3, 'Color', color1, 'LineWidth', 1); hold on;
f1_2 = plot(t, data_com.h(3,:)*1e3, 'Color', color2, 'LineWidth', 1);
f2_1 = plot(t, data_two.h(3,:)*1e3, 'Color', color5, 'LineWidth', 1);
xlabel 'Time (s)'; ylabel '\eta_3 (mJ)';
legend([f1_1,f1_2,f2_1],{'ZohCBF','Log-B','Combined'},'Position',[0.395 0.57 0.23 0.31]);
axis([0 1100 -0.06 0]);
set(gcf, 'Position', [390 300 560 180]);
ax3 = gca;
set(ax3, 'FontSize', 12);
set(ax2, 'Position', ax3.Position);

for i=1:4
figure(3+i); clf;
plot(t, data_nom.u(i,:)*1e3, 'Color', color_u(i,:), 'LineWidth', 1); hold on;
plot(t, data_com.u(i,:)*1e3, '--', 'Color', color_u(i,:), 'LineWidth', 1);
plot(t, data_com.u(i,:)*1e3, ':', 'Color', color_u(i,:), 'LineWidth', 1);
plot([0 1100], [1 1]*wheel_limit*1e3, 'k--');
plot([0 1100], -[1 1]*wheel_limit*1e3, 'k--');
xlabel 'Time (s)'; ylabel(['u_' num2str(i) ' (mNm)']);
legend({'ZohCBF', 'Comparison', 'Combined'}, 'Position', [0.705 0.615 0.20 0.31]);
% axis([0 350 -0.75 0.75]);
set(gcf, 'Position', [1000, 1100-300*(i-1), 650, 200]);
set(gca, 'FontSize', 12);
end

if 0
    figure(1); print -depsc AzElNonCon_legend.eps;
    figure(2); print -depsc ConstraintQNonCon_small.eps;
    figure(3); print -depsc ConstraintENonCon_small.eps;
end