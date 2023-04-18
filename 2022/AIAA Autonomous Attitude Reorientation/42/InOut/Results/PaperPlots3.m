time = load('ZohCBF/time.42');
Safety_nom = load('Nominal/Safety.42');
Safety_zoh = load('ZohCBF/Safety.42');
qbn_nom = load('Nominal/qbn.42'); qbn_nom = qbn_nom(:,[4,1,2,3]);
qbn_zoh = load('ZohCBF/qbn.42'); qbn_zoh = qbn_zoh(:,[4,1,2,3]);
Wheels_nom = load('Nominal/Whl.42');
Wheels_zoh = load('ZohCBF/Whl.42');

color1 = [0.0; 0.5; .3]; % instrument
color2 = [0.85; 0.5; 0]; % tracker
color3 = [0.7; 0.2; 0]; % tracker to moon
color_i = [0.3; 0.95; 0.95]; % instrument sun exclusion color
color_t = [0.95; 0.95; 0]; % tracker sun exclusion color
color_m = [0.9; 0.85; .55]; % tracker moon exlcusion color
color_w = [0.4; 0.4; 0.4];
color_u = [0.000, 0.447, 0.741; 0.850, 0.325, 0.098; 0.929, 0.694, 0.125; 0.494, 0.184, 0.556];    

instrument = [0; 0; 1];
tracker = [sind(10); -cosd(10); 0];
wheel_limit = 0.7e-3;

vectors = [-0.686011 -0.652338 -0.322247;
    -0.994940 -0.099484 -0.014024;
    0.800938 0.368626 -0.471819;
    -0.984159 0.168226 0.055959;
    -0.671108 -0.678372 -0.299041;
    -0.410208 0.836350 0.363658;
    0.997122 -0.039660 -0.064618;
    0.473523 0.809465 0.347193;
    0.999668 0.000141 -0.025778];

svn = [-0.90463166857598032, -0.39102558525948933, -0.16953033911932042]'; % sun vector
mvn = [-0.98433583203767128, 0.16302796725297231, 0.06711819172228313]'; % moon vector

f1 = figure(1); clf;
x1 = RotQ(instrument, qbn_zoh');
x2 = RotQ(tracker, qbn_zoh');
x3 = RotQ(instrument, qbn_nom');
x4 = RotQ(tracker, qbn_nom');
az1 = atan2d(x1(2,:), x1(1,:)); az1(az1 < 0) = az1(az1<0) + 360;
az2 = atan2d(x2(2,:), x2(1,:)); az2(az2 < 0) = az2(az2<0) + 360;
az3 = atan2d(x3(2,:), x3(1,:)); az3(az3 < 0) = az3(az3<0) + 360;
az4 = atan2d(x4(2,:), x4(1,:)); az4(az4 < 0) = az4(az4<0) + 360;
az1(logical([0, abs(diff(az1)) > 100])) = NaN;
az2(logical([0, abs(diff(az2)) > 100])) = NaN;
az3(logical([0, abs(diff(az3)) > 100])) = NaN;
az4(logical([0, abs(diff(az4)) > 100])) = NaN;
el1 = asind(x1(3,:));
el2 = asind(x2(3,:));
el3 = asind(x3(3,:));
el4 = asind(x4(3,:));

region_axes = axes(f1, 'FontSize', 14);
[c1, c2, c3] = PlotCone(svn, deg2rad(45));
indices = abs([c1(:), c2(:), c3(:)]*svn - cosd(45)) < 1e-5;
azc = atan2d(c2(indices), c1(indices)); azc(azc < 0) = azc(azc<0) + 360;
elc = asind(c3(indices));
region1 = fill(azc,elc,color_t','EdgeAlpha',0); hold on;

[c1, c2, c3] = PlotCone(mvn, deg2rad(30));
indices = abs([c1(:), c2(:), c3(:)]*mvn - cosd(30)) < 1e-5;
azc = atan2d(c2(indices), c1(indices));
azc(azc < 0) = azc(azc<0) + 360;
elc = asind(c3(indices));
region2 = fill(azc,elc,color_m','EdgeAlpha',0); hold on;

[c1, c2, c3] = PlotCone(svn, deg2rad(25));
indices = abs([c1(:), c2(:), c3(:)]*svn - cosd(25)) < 1e-5;
azc = atan2d(c2(indices), c1(indices));
azc(azc < 0) = azc(azc<0) + 360;
elc = asind(c3(indices));
region3 = fill(azc,elc,color_i','EdgeAlpha',0); hold on;

line_axes = copyobj(region_axes, f1);
delete(get(line_axes, 'Children'));
f1 = plot(az1, el1, 'LineWidth', 3, 'Color', color1); hold on;
f2 = plot(az2, el2, 'LineWidth', 3, 'Color', color2);
f3 = plot(az3, el3, '--', 'LineWidth', 3, 'Color', color1);
f4 = plot(az4, el4, '--', 'LineWidth', 3, 'Color', color2);
azs = atan2d(vectors(:,2), vectors(:,1)); azs(azs < 0) = azs(azs<0) + 360;
els = asind(vectors(:,3));
f5 = plot(azs, els, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(az1(1), el1(1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
plot(az2(1), el2(1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
xlabel('Azimuth in Inertial Frame (deg)','FontSize',15); ylabel('Elevation in Inertial Frame (deg)', 'FontSize', 15);
legend(line_axes,[f1,f2,f3,f4],{'ZohCBF b_1','ZohCBF b_2','PD b_1','PD b_2'},...
    'FontSize',11,'Position',[0.593 0.669 0.20 0.264]); %,'Color','none');
legend(region_axes,[region3,region1,region2,f5],{'b_1 Sun Keep-out','b_2 Sun Keep-out','b_2 Moon Keep-out','Targets for b_1'},...
    'FontSize',11,'Position',[0.131 0.634 0.229 0.249],'Color','none');
set(gcf, 'Position', [200 400 750 400]);

axis(region_axes, 'equal');
axis(region_axes, [0, 360, -90, 90]);
axis(line_axes, 'equal');
axis(line_axes, [0, 360, -90, 90]);
set(line_axes, 'Color', 'none');

figure(2); clf;
mu = 0.00167;
h_zoh = Safety_zoh(:,[1,3,5]);
hdot_zoh = Safety_zoh(:,[2,4,6]);
H_zoh = h_zoh + absSq(hdot_zoh)/(2*mu);
h_nom = Safety_nom(:,[1,3,5]);
hdot_nom = Safety_nom(:,[2,4,6]);
H_nom = h_nom + absSq(hdot_nom)/(2*mu);
f1 = plot(time, H_zoh(:,1), 'Color', color1, 'LineWidth', 1); hold on;
f2 = plot(time, H_zoh(:,2), 'Color', color2, 'LineWidth', 1);
f3 = plot(time, H_zoh(:,3), 'Color', color3, 'LineWidth', 1);
f4 = plot(time, H_nom(:,1), '--', 'Color', color1, 'LineWidth', 1); hold on;
f5 = plot(time, H_nom(:,2), '--', 'Color', color2, 'LineWidth', 1);
f6 = plot(time, H_nom(:,3), '--', 'Color', color3, 'LineWidth', 1);
xlabel 'Time (s)'; ylabel 'h_{ }';
legend([f1,f2,f3],{'h_{instr-sun}','h_{track-sun}','h_{track-moon}'},...
    'Orientation','horizontal','Position',[0.322,0.245,0.55,0.1475]);
axis([0 2100 -2.7 0.4]);
set(gcf, 'Position', [390 550 560 220]);
ax2 = gca;
set(ax2, 'FontSize', 12);

figure(3); clf;
f1 = plot(time, Safety_zoh(:,7)*1e3, 'Color', color_w, 'LineWidth', 1); hold on;
f2 = plot(time, Safety_nom(:,7)*1e3, '--', 'Color', color_w, 'LineWidth', 1);
xlabel 'Time (s)'; ylabel '\eta_4 (mJ)';
legend([f1,f2],{'ZohCBF','PD'},'Position',[0.725,0.27,0.16,0.21]);
axis([0 2100 -0.2 0.18]);
set(gcf, 'Position', [390 250 560 200]);
ax3 = gca;
set(ax3, 'FontSize', 12);
set(ax2, 'Position', [ax3.Position(1), ax2.Position(2), ax3.Position(3), ax2.Position(4)]);

for i=1:4
figure(3+i); clf;
plot(time, Wheels_zoh(:,(i-1)*3+2)*1e3, 'Color', color_u(i,:), 'LineWidth', 1); hold on;
plot(time, Wheels_nom(:,(i-1)*3+2)*1e3, '--', 'Color', color_u(i,:), 'LineWidth', 1);
plot([0 2100], [1 1]*wheel_limit*1e3, 'k--');
plot([0 2100], -[1 1]*wheel_limit*1e3, 'k--');
xlabel 'Time (s)'; ylabel(['u_' num2str(i) ' (mNm)']);
legend ZohCBF PD Location North;
axis([0 2100 -0.75 0.75]);
set(gcf, 'Position', [1000, 1100-300*(i-1), 650, 200]);
set(gca, 'FontSize', 12);
end

if 0
    figure(1); print -depsc AzEl_42_legend.eps;
    figure(2); print -depsc ConstraintQ_42_small.eps;
    figure(3); print -depsc ConstraintE_42.eps;
    figure(4); print -depsc Control1_42.eps;
    figure(5); print -depsc Control2_42.eps;
    figure(6); print -depsc Control3_42.eps;
    figure(7); print -depsc Control4_42.eps;
end