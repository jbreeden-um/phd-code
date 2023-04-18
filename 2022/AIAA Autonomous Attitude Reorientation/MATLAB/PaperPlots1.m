%% Regular Cases
data1 = load('Results/CBF.mat');
data2 = load('Results/Barrier.mat');
data3 = load('Results/SMC.mat');
data4 = load('Results/MPC.mat');

try
    cross(get_s(0),[1;0;0]);
catch
    ComputeConstants;
end

color1 = [0.2; 0.2; 1];
color2 = [0; 0.6; 0];
color3 = [0.7; 0; 0];
color4 = [0.5; 0.5; 0.5];
color_w = [0.4; 0.4; 0.4];
color_u = [0.000, 0.447, 0.741; 0.850, 0.325, 0.098; 0.929, 0.694, 0.125; 0.494, 0.184, 0.556];    

t = data1.t;
wheel_limit = data1.constants.wheel_limit;
p1 = data1.constants.p1;
p2 = data1.constants.p2;
ctheta = data1.ctheta;
s_target = data1.s_target;

%%
figure(1); clf;
x1_1 = RotQ(p1, data1.x(:,1:4)');
x1_2 = RotQ(p2, data1.x(:,1:4)');
x2_1 = RotQ(p1, data2.x(:,1:4)');
x2_2 = RotQ(p2, data2.x(:,1:4)');
x3_1 = RotQ(p1, data3.x(:,1:4)');
x3_2 = RotQ(p2, data3.x(:,1:4)');
x4_1 = RotQ(p1, data4.x(:,1:4)');
x4_2 = RotQ(p2, data4.x(:,1:4)');
az1_1 = atan2d(x1_1(2,:), x1_1(1,:));
az1_2 = atan2d(x1_2(2,:), x1_2(1,:));
az2_1 = atan2d(x2_1(2,:), x2_1(1,:));
az2_2 = atan2d(x2_2(2,:), x2_2(1,:));
az3_1 = atan2d(x3_1(2,:), x3_1(1,:));
az3_2 = atan2d(x3_2(2,:), x3_2(1,:));
az4_1 = atan2d(x4_1(2,:), x4_1(1,:));
az4_2 = atan2d(x4_2(2,:), x4_2(1,:));
az1_2(logical([0, abs(diff(az1_2)) > 1])) = NaN; % to prevent discontinuities
az2_2(logical([0, abs(diff(az2_2)) > 1])) = NaN;
az3_2(logical([0, abs(diff(az3_2)) > 1])) = NaN;
az4_2(logical([0, abs(diff(az4_2)) > 1])) = NaN;
el1_1 = asind(x1_1(3,:));
el1_2 = asind(x1_2(3,:));
el2_1 = asind(x2_1(3,:));
el2_2 = asind(x2_2(3,:));
el3_1 = asind(x3_1(3,:));
el3_2 = asind(x3_2(3,:));
el4_1 = asind(x4_1(3,:));
el4_2 = asind(x4_2(3,:));
[c1, c2, c3] = PlotCone(get_s(0), acos(ctheta));
indices = abs([c1(:), c2(:), c3(:)]*get_s(0) - ctheta) < 1e-5;
azc = atan2d(c2(indices), c1(indices));
elc = asind(c3(indices));
fill(azc,elc,'r','EdgeAlpha',0,'FaceColor',[1;.2;.2]); hold on;
f1_1 = plot(az1_1, el1_1, '-', 'LineWidth', 3, 'Color', color1);
f1_2 = plot(az1_2, el1_2, '--', 'LineWidth', 3, 'Color', color1);
f2_1 = plot(az2_1, el2_1, '-', 'LineWidth', 3, 'Color', color2);
f2_2 = plot(az2_2, el2_2, '--', 'LineWidth', 3, 'Color', color2);
f3_1 = plot(az3_1, el3_1, '-', 'LineWidth', 3, 'Color', color3);
f3_2 = plot(az3_2, el3_2, '--', 'LineWidth', 3, 'Color', color3);
f4_1 = plot(az4_1, el4_1, '-', 'LineWidth', 3, 'Color', color4);
f4_2 = plot(az4_2, el4_2, '--', 'LineWidth', 3, 'Color', color4);
azs = atan2d(s_target(2), s_target(1));
els = asind(s_target(3));
f5 = plot(azs, els, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
f6 = plot(az1_1(1), el1_1(1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
f7 = plot(az1_2(1), el1_2(1), 'ko', 'LineWidth', 3, 'MarkerSize', 8);
xlabel 'Azimuth in Inertial Frame (deg)'; ylabel 'Elevation in Inertial Frame (deg');
l = legend([f1_1,f2_1,f3_1,f4_1,f6,f7,f5],{'ZohCBF','Log-B','SMC','NMPC','Initial b_1','Initial b_2','Target for b_1'},'FontSize',12);
axis equal;
axis([-180, 180, -90, 90]);
set(gcf, 'Position', [200 900 750 400]);
set(gca, 'FontSize', 14);
set(l, 'Position', [0.695555559171571 0.444 0.189333329717318 0.437499987483025])

%%
figure(2); clf;
f1_1 = plot(data1.t, data1.h(1,:), '-', 'Color', color1, 'LineWidth', 1); hold on;
f1_2 = plot(data1.t, data1.h(2,:), '--', 'Color', color1, 'LineWidth', 1);
f2_1 = plot(data2.t, data2.h(1,:), '-', 'Color', color2, 'LineWidth', 1);
f2_2 = plot(data2.t, data2.h(2,:), '--', 'Color', color2, 'LineWidth', 1);
f3_1 = plot(data3.t, data3.h(1,:), '-', 'Color', color3, 'LineWidth', 1);
f3_2 = plot(data3.t, data3.h(2,:), '--', 'Color', color3, 'LineWidth', 1);
f4_1 = plot(data4.t, data4.h(1,:), '-', 'Color', color4, 'LineWidth', 1);
f4_2 = plot(data4.t, data4.h(2,:), '--', 'Color', color4, 'LineWidth', 1);
f1 = plot(0,0,'k-', 'LineWidth', 1);
f2 = plot(0,0,'k--', 'LineWidth', 1);
xlabel 'Time (s)'; ylabel '\kappa_{ }';
% legend([f1,f2,f3,f4],{'CBF H_1','CBF H_2','L.F. H_1','L.F. H_2'},...
%     'Orientation','horizontal','Position',[0.145 0.29 0.75 0.15])
l = legend([f1,f2],{'\kappa_1','\kappa_2'},'Orientation','horizontal','Location','SouthEast');
axis([0 500 -1.4 0.1]);
set(gcf, 'Position', [390 600 560 180]); %, 'Renderer', 'Painters');
ax2 = gca;
set(ax2, 'FontSize', 12);
set(l, 'Position', [0.494047620750609 0.651851862183323 0.233928569725581 0.163888884749678]);

figure(3); clf;
f1 = plot(data1.t, data1.h(3,:)*1e3, 'Color', color1, 'LineWidth', 1); hold on;
f2 = plot(data2.t, data2.h(3,:)*1e3, 'Color', color2, 'LineWidth', 1);
f3 = plot(data3.t, data3.h(3,:)*1e3, 'Color', color3, 'LineWidth', 1);
f4 = plot(data4.t, data4.h(3,:)*1e3, 'Color', color4, 'LineWidth', 1);
xlabel 'Time (s)'; ylabel '\eta_3 (mJ)';
l = legend([f1,f2,f3,f4],{'ZohCBF','Log-B','SMC','NMPC'});
axis([0 500 -0.06 0]);
set(gcf, 'Position', [390 300 560 180]); %, 'Renderer', 'Painters');
ax3 = gca;
set(ax3, 'FontSize', 12);
set(ax2, 'Position', ax3.Position);
set(l, 'Position', [0.685119051450774 0.473148160731350 0.219642853311130 0.452777765194575]);

%%
for i=1:4
figure(3+i); clf;
plot(data1.t, data1.u(i,:)*1e3, 'Color', color1, 'LineWidth', 1); hold on;
plot(data2.t, data2.u(i,:)*1e3, 'Color', color2, 'LineWidth', 1);
plot(data3.t, data3.u(i,:)*1e3, 'Color', color3, 'LineWidth', 1);
plot(data4.t, data4.u(i,:)*1e3, 'Color', color4, 'LineWidth', 1);
plot([0 500], [1 1]*wheel_limit*1e3, 'k--', 'LineWidth', 1);
plot([0 500], -[1 1]*wheel_limit*1e3, 'k--', 'LineWidth', 1);
xlabel 'Time (s)'; ylabel(['u_' num2str(i) ' (mNm)']);
legend ZohCBF 'Log-B' SMC NMPC Location East;
axis([0 500 -0.75 0.75]);
set(gcf, 'Position', [1000, 1100-300*(i-1), 650, 170]);
set(gca, 'FontSize', 12);
end

if 0
    figure(1); print -depsc AzEl_legend.eps;
    figure(2); print -depsc ConstraintQ_small.eps;
    figure(3); print -depsc ConstraintE_small.eps;
    figure(4); print -depsc Control1.eps;
    figure(5); print -depsc Control2.eps;
    figure(6); print -depsc Control3.eps;
    figure(7); print -depsc Control4.eps;
end