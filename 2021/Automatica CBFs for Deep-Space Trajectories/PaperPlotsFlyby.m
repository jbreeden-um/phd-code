data1 = load('Results/SimConstant.txt');
data2 = load('Results/SimVariable.txt');
data3 = load('Results/SimIntegrated2S.txt');
data4 = load('Results/SimIntegrated3S.txt');

control1 = load('Results/ControlConstant.txt');
control2 = load('Results/ControlVariable.txt');
control3 = load('Results/ControlIntegrated2S.txt');
control4 = load('Results/ControlIntegrated3S.txt');

t1 = data1(:,1);      t2 = data2(:,1);      t3 = data3(:,1);      t4 = data4(:,1);
r1 = data1(:,2:4);    r2 = data2(:,2:4);    r3 = data3(:,2:4);    r4 = data4(:,2:4);
v1 = data1(:,5:7);    v2 = data2(:,5:7);    v3 = data3(:,5:7);    v4 = data4(:,5:7);
% H1 = data1(:,8);      H2 = data2(:,8);      H3 = data3(:,8);      H4 = data4(:,8);
H1 = control1(:,2);   H2 = control2(:,2);   H3 = control3(:,2);   H4 = control4(:,2);
% u1 = data1(:,9:11);   u2 = data2(:,9:11);   u3 = data3(:,9:11);   u4 = data4(:,9:11);
u1 = control1(:,3:5); u2 = control2(:,3:5); u3 = control3(:,3:5); u4 = control4(:,3:5);
sigma1 = data1(:,12); sigma2 = data2(:,12); sigma3 = data3(:,12); sigma4 = data4(:,12);
                                            beta3 = data3(:,13);  beta4 = data3(:,13);

color1 = [0;0;1];
color2 = [0;1;0];
color3 = [1;0;0];
color4 = [1;0;1];

day = 24*3600;
t1 = t1/day; t2 = t2/day; t3 = t3/day; t4 = t4/day;
te = 6e6/day;
u1 = u1*1e3; u2 = u2*1e3; u3 = u3*1e3; u4 = u4*1e3;

%%
figure(1); clf;
[s1, s2, s3] = sphere(20);
r_Ceres = 476000; 
hold on;
p1 = plot3(r1(:,1)/1e6, r1(:,2)/1e6, r1(:,3)/1e6, 'LineWidth', 2, 'Color', color1);
p2 = plot3(r2(:,1)/1e6, r2(:,2)/1e6, r2(:,3)/1e6, 'LineWidth', 2, 'Color', color2);
p3 = plot3(r3(:,1)/1e6, r3(:,2)/1e6, r3(:,3)/1e6, 'LineWidth', 2, 'Color', color3);
p4 = plot3(r4(:,1)/1e6, r4(:,2)/1e6, r4(:,3)/1e6, 'LineWidth', 2, 'Color', color4);
xlabel 'X (Mm)';
ylabel 'Y (Mm)';
zlabel 'Z (Mm)';
plot3(0, 0, 0, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', [0.5;0.5;0.5]);
axis equal;
set(gcf, 'Position', [1700 800 600 500])
axis([-5.5e7 5.5e7 -5e7 5e7 -2e7 2e7]/1e6);
view([0;0;1])
legend([p1,p2,p3,p4],{'Case $H^A_1$','Case $H^A_2$','Case $H^A_3$','Case $H^A_4$'},'Location','East','Interpreter','latex','FontSize',14)
grid on;

%%
surf(s1*r_Ceres/1e6, s2*r_Ceres/1e6, s3*r_Ceres/1e6, 'FaceColor', [0.7;0.7;0.7], 'FaceAlpha', 1);
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
set(gca,'FontSize',12)
legend([p4],{'Case $H^A_4$'},'Location','SouthEast','Interpreter','latex','FontSize',18)

%%
figure(2); clf;
p1 = semilogy(t1, H1, 'Color', color1); hold on;
p2 = semilogy(t2, H2, 'Color', color2);
p3 = semilogy(t3, H3, 'Color', color3);
p4 = semilogy(t4, H4, 'Color', color4);
xlabel 'Time (days)';
ylabel 'H (m)';
legend([p1,p2,p3,p4],{'Case $H_1$', 'Case $H_2$', 'Case $H_3$', 'Case $H_4$'},...
    'Location','NorthEast','Interpreter','latex','FontSize',13)
set(gcf, 'Position', [1850 500 560 200]);
axis([0 te -1e9 -1e3])
set(gca, 'YTick', [-1e8 -1e7 -1e6 -1e5 -1e4], 'FontSize', 11)

%%
figure(3); clf;
p1 = plot(t1,u1(:,1),'Color',color1); hold on;
p2 = plot(t2,u2(:,1),'Color',color2);
xlabel 'Time (days)';
ylabel 'u_x (mm/s^2)';
hold on;
plot([t1(1), t1(end)], [1e-1 1e-1], 'k--');
plot([t1(1), t1(end)], [-1e-1 -1e-1], 'k--');
set(gca, 'FontSize', 11);
a = gca;
legend([p1,p2],{'Case $H_1$','Case $H_2$'}, 'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', 13)
b = copyobj(a,gcf);
delete(get(b,'Children'));
set(b,'Color','none','XTick',0, 'YTick',[])%,'XTick',[],'YTick',[],'Box','off')
set(b.XLabel,'String', ' ','Color',[0;0;0]);
set(b.YLabel,'String', ' ');
set(b, 'XTickLabel', {' '});
p3 = plot(t3,u3(:,1),'Color',color3,'Parent',b);
p4 = plot(t4,u4(:,1),'Color',color4,'Parent',b);
legend([p3,p4],{'Case $H_3$','Case $H_4$'}, 'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 13, 'Color', 'w')
set(gcf, 'Position', [2550 800 560 180]);
set(a.Legend, 'Position', [0.70 0.64 0.20 0.26]);
set(b.Legend, 'Position', [0.70 0.29 0.20 0.26]);
axis(a, [0 te -1.05e-1 1.05e-1])
axis(b, [0 te -1.05e-1 1.05e-1])

figure(4); clf;
p1 = plot(t1,u1(:,2),'Color',color1); hold on;
p2 = plot(t2,u2(:,2),'Color',color2);
% p3 = plot(t3,u3(:,2),'Color',color3);
% p4 = plot(t4,u4(:,2),'Color',color4);
xlabel 'Time (days)';
ylabel 'u_y (mm/s^2)';
hold on;
plot([t1(1), t1(end)], [1e-1 1e-1], 'k--');
plot([t1(1), t1(end)], [-1e-1 -1e-1], 'k--');
set(gca, 'FontSize', 11);
a = gca;
legend([p1,p2],{'Case $H_1$','Case $H_2$'}, 'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', 13)
b = copyobj(a,gcf);
delete(get(b,'Children'));
set(b,'Color','none','XTick',0, 'YTick',[])%,'XTick',[],'YTick',[],'Box','off')
set(b.XLabel,'String', ' ','Color',[0;0;0]);
set(b.YLabel,'String', ' ');
set(b, 'XTickLabel', {' '});
p3 = plot(t3,u3(:,2),'Color',color3,'Parent',b);
p4 = plot(t4,u4(:,2),'Color',color4,'Parent',b);
legend([p3,p4],{'Case $H_3$','Case $H_4$'}, 'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 13, 'Color', 'w')
set(gcf, 'Position', [2550 500 560 180]);
set(a.Legend, 'Position', [0.70 0.64 0.20 0.26]);
set(b.Legend, 'Position', [0.70 0.29 0.20 0.26]);
axis(a, [0 te -1.05e-1 1.05e-1])
axis(b, [0 te -1.05e-1 1.05e-1])

figure(5); clf;
p1 = plot(t1,u1(:,3),'Color',color1); hold on;
p2 = plot(t2,u2(:,3),'Color',color2);
% p3 = plot(t3,u3(:,3),'Color',color3);
% p4 = plot(t4,u4(:,3),'Color',color4);
xlabel 'Time (days)';
ylabel 'u_z (mm/s^2)';
hold on;
plot([t1(1), t1(end)], [1e-1 1e-1], 'k--');
plot([t1(1), t1(end)], [-1e-1 -1e-1], 'k--');
set(gca, 'FontSize', 11);
a = gca;
legend([p1,p2],{'Case $H_1$','Case $H_2$'}, 'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', 13)
b = copyobj(a,gcf);
delete(get(b,'Children'));
set(b,'Color','none','XTick',0, 'YTick',[])%,'XTick',[],'YTick',[],'Box','off')
set(b.XLabel,'String', ' ','Color',[0;0;0]);
set(b.YLabel,'String', ' ');
set(b, 'XTickLabel', {' '});
p3 = plot(t3,u3(:,3),'Color',color3,'Parent',b);
p4 = plot(t4,u4(:,3),'Color',color4,'Parent',b);
legend([p3,p4],{'Case $H_3$','Case $H_4$'}, 'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 13, 'Color', 'w')
set(gcf, 'Position', [2550 200 560 180]);
set(a.Legend, 'Position', [0.70 0.64 0.20 0.26]);
set(b.Legend, 'Position', [0.70 0.29 0.20 0.26]);
axis(a, [0 te -1.05e-1 1.05e-1])
axis(b, [0 te -1.05e-1 1.05e-1])

%%
figure(6); clf;
p1 = semilogy(t1, -r_Ceres + vecnorm(r1,2,2), 'Color', color1); hold on;
p2 = semilogy(t2, -r_Ceres + vecnorm(r2,2,2), 'Color', color2);
p3 = semilogy(t3, -r_Ceres + vecnorm(r3,2,2), 'Color', color3);
p4 = semilogy(t4, -r_Ceres + vecnorm(r4,2,2), 'Color', color4);
xlabel 'Time (days)';
ylabel 'Altitude Above Ceres (m)';
legend([p1,p2,p3,p4],{'Case $H_1$', 'Case $H_2$', 'Case $H_3$', 'Case $H_4$'},...
    'Location', 'SouthEast' ,'Interpreter', 'latex', 'FontSize', 13)
set(gcf, 'Position', [1850 200 560 200]);
axis([0 te 1e5 1e9]);
set(gca, 'YTick', [1e5 1e6 1e7 1e8 1e9], 'FontSize', 11)

%%
data5 = load('Results/SimConstantNoSwitch.txt');
t5 = data5(:,1)/day; H5 = data5(:,8);
color5 = [0;0;0];
figure(7); clf;
semilogy(t1, H1, 'Color', color2); hold on;
semilogy(t5, H5, 'Color', color5);
xlabel 'Time (days)';
ylabel 'H (m)';
legend 'Switched' 'Always Active' Location NorthEast Interpreter latex
set(gcf, 'Position', [1850 800 560 200]);
axis([0 te -1e9 -1e-4])
set(gca, 'FontSize', 11)
l = legend;
set(l, 'FontSize', 13, 'Interpreter', 'latex');
yticks([-1e5,-1])

%%
figure(8); clf;
tpts = [0, 2, 3, 5, 6.5, 10];
Hpts = [-3, -1, -0.1, -1, -2, -1.5];
ts = linspace(tpts(1), tpts(end), 1000);
Hs = spline(tpts, Hpts, ts);
sigma = Hs > -1 | (Hs > -2 & ts>3 & ts<8);
figure(8); clf;
subplot(2,1,1); plot(ts, sigma, 'b', 'LineWidth', 3); ylabel 'Discrete State (\sigma)'
subplot(2,1,2); 
plot(ts, Hs, 'b', 'LineWidth', 2); hold on;
plot(ts, -ones(size(ts)), 'r--');
plot(ts, -2*ones(size(ts)), 'g--');
xlabel 'Time'; ylabel 'RCBF (H)';

[~,i1] = min(abs(Hs(1:300) + 2));
[~,i2] = min(abs(Hs(1:300) + 1));
[~,i3] = min(abs(Hs(301:800) + 1)); i3 = i3 + 300;
[~,i4] = min(abs(Hs(301:800) + 2)); i4 = i4 + 300;
[~,i5] = min(abs(Hs(801:end) + 2)); i5 = i5 + 800;
plot([ts(i1), ts(i1)], [-3 0], 'k--');
plot([ts(i2), ts(i2)], [-3 0], 'k--');
plot([ts(i3), ts(i3)], [-3 0], 'k--');
plot([ts(i4), ts(i4)], [-3 0], 'k--');
plot([ts(i5), ts(i5)], [-3 0], 'k--');
legend({'$H$','$-\epsilon_1$','$-\epsilon_2$'},'Interpreter','latex',...
    'FontSize', 12, 'Position',[0.696 0.275 0.125 0.15]);

subplot(2,1,1); hold on;
plot([ts(i1), ts(i1)], [1 0], 'k--');
plot([ts(i2), ts(i2)], [1 0], 'k--');
plot([ts(i3), ts(i3)], [1 0], 'k--');
plot([ts(i4), ts(i4)], [1 0], 'k--');
plot([ts(i5), ts(i5)], [1 0], 'k--');
set(gcf, 'Position', [1910 845 560 300])