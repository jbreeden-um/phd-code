data = load('Results/SimEros.txt');
control = load('Results/ControlEros.txt');

t = data(:,1);
r = data(:,2:4);
v = data(:,5:7);
% H = data(:,8);
H = control(:,2);
% u = data(:,9:11);
u = control(:,3:5);
sigma = data(:,12);
te = data(:,13);

rho = 500;
epsilon1 = 100;
figure(1); clf;
plot(t, H);
xlabel 'Time (s)'; ylabel 'Max CBF Value (m)';
set(gcf, 'Position', [1850 950 560 200]);
set(gca, 'FontSize', 11);

Eros = load('Eros.mat');
r_rot = zeros(size(r));
omega = [0.0003101010326521003; 6.231917786350903e-5; 9.81048895288355e-5];
for i=1:length(t)
    rot = expm(skew(t(i)*omega)); % body to inertial
    s = rot'*r(i,1:3)';
    r_rot(i,:) = s'/1e3;
end
figure(2); clf;
trisurf(Eros.plates+1, Eros.vertices(:,1), Eros.vertices(:,2), Eros.vertices(:,3), ...
    'FaceColor', [.7; .7; .7], 'FaceAlpha', 1, 'EdgeAlpha', 0.3);
hold on;
plot3(r_rot(:,1), r_rot(:,2), r_rot(:,3),'LineWidth',4); hold on;
plot3(r_rot(1,1), r_rot(1,2), r_rot(1,3),'ko','MarkerFaceColor','b');
r_final = rot'*[20;0;0];
plot3(r_final(1),r_final(2),r_final(3),'ko','MarkerFaceColor','g');
xlabel 'r_x (km)';
ylabel 'r_y (km)';
zlabel 'r_z (km)';
axis equal;
view(-50,30);
set(gcf, 'Position', [1800 100 800 600])
set(gca, 'FontSize', 11); 

figure(3); clf;
plot(t,u(:,1));
xlabel 'Time (s)';
ylabel 'u_x (m/s^2)';
hold on;
plot([t(1), t(end)], [0.1 0.1], 'r--');
plot([t(1), t(end)], [-0.1 -0.1], 'r--');
axis([0 6000 -0.105 0.105])
set(gcf, 'Position', [2650 900 560 160]);
set(gca, 'FontSize', 11); 

figure(4); clf;
plot(t,u(:,2));
xlabel 'Time (s)';
ylabel 'u_y (m/s^2)';
hold on;
plot([t(1), t(end)], [0.1 0.1], 'r--');
plot([t(1), t(end)], [-0.1 -0.1], 'r--');
axis([0 6000 -0.105 0.105])
set(gcf, 'Position', [2650 600 560 160]);
set(gca, 'FontSize', 11); 

figure(5); clf;
plot(t,u(:,3));
xlabel 'Time (s)';
ylabel 'u_z (m/s^2)';
hold on;
plot([t(1), t(end)], [0.1 0.1], 'r--');
plot([t(1), t(end)], [-0.1 -0.1], 'r--');
axis([0 6000 -0.105 0.105])
set(gcf, 'Position', [2650 300 560 160]);
set(gca, 'FontSize', 11); 

return
%% Intro Figure
fig = figure(6); clf;
plotEros = trisurf(Eros.plates+1, Eros.vertices(:,1), Eros.vertices(:,2), Eros.vertices(:,3), ...
    'FaceColor', [.5; .5; .5], 'FaceAlpha', 1, 'EdgeAlpha', 0.3);
hold on;
r_final = [20e3;0;0];
xlabel 'r_x (km)'; ylabel 'r_y (km)'; zlabel 'r_z (km)';
axis equal;
set(gcf, 'Position', [1500 100 1920 1080])
set(gcf, 'color', 'k')
axis off
set(gcf, 'InvertHardCopy', 'off');

n = 200;
points = zeros(3,n);
theta = linspace(-pi, 0, n);
v = Eros.vertices;
polar = [vecnorm(v,2,2), atan2(v(:,2), v(:,1)), asin(v(:,3)./vecnorm(v(:,1:2),2,2))];
polar = polar(abs(polar(:,3)) < .1, :);
for i=1:n
    vals = v(:,1:3)*[cos(theta(i)); sin(theta(i)); -0.2*sin(theta(i))]./vecnorm(v(:,1:3),2,2);
    [~,index] = max(vals);
    radius = norm(v(index,:))+2+sin(theta(i));
    points(:,i) = radius*[cos(theta(i)); sin(theta(i)); 0];
end
p1 = polyfit(theta, points(1,:), 8);

[xData, yData] = prepareCurveData(theta, points(2,:));
ft = fittype( 'sin6' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf];
opts.StartPoint = [7.77158051636564 1 -0.0178081403573535 2.89447014862312 2 -1.08194602093731 1.01595765051632 4 1.286724202919 0.716020575148793 6 2.29027389692948 0.588844132982594 8 2.52230600376338 0.197865342200652 10 2.5379272593807];
[fity, ~] = fit( xData, yData, ft, opts );

theta_plot = linspace(-3, -0.53, 1000);
x = polyval(p1, theta_plot);
y = fity(theta_plot);
z = 0*theta_plot;
plot3(x, y, z, 'Color', [0 0.447 0.741], 'LineWidth', 6);
view(-5,90)
plot3(x(1), y(1), z(1),'wo','MarkerFaceColor',[0 0.447 0.741],'MarkerSize',10,'LineWidth',2);
plot3(x(end), y(end), z(end),'wo','MarkerFaceColor',[0 0.447 0.741],'MarkerSize',10,'LineWidth',2);

%%%
try
    delete(v_line)
    delete(v_arrow)
catch
end
theta_sc = -2.36;
xdot = polyval(polyder(p1), theta_sc);
[ydot,~] = differentiate(fity, theta_sc);
zdot = 0;
vector = [xdot; ydot; zdot]; vector = vector/norm(vector);
R = get_rotation([0;0;1], vector);
len = 4;
radius = 0.3;
height = 0.6;
x0 = polyval(p1, theta_sc);
y0 = fity(theta_sc);
z0 = 0;
[ax,ay,az] = cylinder(1,40);
[cx, cy, cz] = rotate_coordinates(R, ax*radius.*[0; 1], ay*radius.*[0;1], top-0.1-az*height);
v_line = plot3([x0, x0+vector(1)*len], [y0, y0+vector(2)*len], [z0, z0+vector(3)*len], ...
    'Color', [0; 0.6; 0], 'LineWidth', 5);
v_arrow = surf(x0+vector(1)*len+cx, y0+vector(2)*len+cy, z0+vector(3)*len+cz, ...
    'FaceColor', [0; 0.6; 0], 'EdgeAlpha', 0, 'FaceAlpha', 1);

%%%
try
    delete(r_line)
    delete(r_arrow)
catch
end
plot3(0, 0, 11,'wo','MarkerFaceColor',[0 0 0],'MarkerSize',12,'LineWidth',0.5);
r_line = plot3([0, x0+0.15], [0, y0+0.25], [10, z0], 'Color', [1; 0.4; 0], 'LineWidth', 5);
vector = [x0+0.15; y0+0.25; 0];
R = get_rotation([0;0;1], vector);
[cx, cy, cz] = rotate_coordinates(R, ax*radius.*[0; 1], ay*radius.*[0;1], top-0.1-az*height);
r_arrow = surf(vector(1)+cx, vector(2)+cy, vector(3)+cz, 'FaceColor', [1; 0.4; 0], 'EdgeAlpha', 0, 'FaceAlpha', 1);
frame = gca;

%%%
try
    delete(im);
catch
end
sc = axes('pos',[0.446, 0.318, 0.06 0.06]);
img = imread('NEAR.png');
im = imshow(img);
alpha = (ones(size(img,1),size(img,2)));
alpha(sum(img,3)==0) = 0;
set(im, 'AlphaData', alpha);

%%%
try
    delete(t1);
    delete(t2);
catch
end
axes(frame)
t1 = text(-3.2, -1.8, 8, '$r$', 'Interpreter', 'latex', 'Color', [1; 0.6; 0], 'FontSize', 60);
t2 = text(-3.1, -6.0, 8, '$v$', 'Interpreter', 'latex', 'Color', [0; 0.8; 0], 'FontSize', 60);
axes(sc)

return
%% Video
fig = figure(2); clf;
plotEros = trisurf(Eros.plates+1, Eros.vertices(:,1), Eros.vertices(:,2), Eros.vertices(:,3), ...
    'FaceColor', [.5; .5; .5], 'FaceAlpha', 1, 'EdgeAlpha', 0.3);
hold on;
plot3(r(1,1)/1e3, r(1,2)/1e3, r(1,3)/1e3,'ko','MarkerFaceColor','r','MarkerSize',10);
r_final = [20e3;0;0];
plot3(r_final(1)/1e3,r_final(2)/1e3,r_final(3)/1e3,'ko','MarkerFaceColor','g','MarkerSize',10);
xlabel 'r_x (km)'; ylabel 'r_y (km)'; zlabel 'r_z (km)';
axis equal; view(-70,30);
set(gcf, 'Position', [1500 100 1920 1080])
set(gcf, 'color', 'k')
axis off
plotSC = plot3(r(1,1)/1e3, r(1,2)/1e3, r(1,3)/1e3,'ko','MarkerFaceColor',[.6;.6;1],'MarkerSize',10);
set(gcf, 'InvertHardCopy', 'off');

switches = logical(load('Results/SwitchingEros.txt'));
% splot = plot3(NaN, NaN, NaN, 'w.', 'MarkerSize', 35);
splot = plot3(NaN, NaN, NaN, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'w');

record = 1;
video_duration = 60;
frame_rate = 30;
nframes = video_duration*frame_rate;
moviename = 'MovieEros';
vidObj = VideoWriter(moviename);
vidObj.FrameRate = frame_rate;
if record, open(vidObj); end

for i=1:nframes
    tframe = i*t(end)/nframes;
    rframe = interp1(t, r, tframe)/1e3;
    rot = expm(skew(tframe*omega)); % body to inertial
    vertices = Eros.vertices*rot';
    set(plotEros, 'Vertices', vertices);
    set(plotSC, 'XData', rframe(1), 'YData', rframe(2), 'ZData', rframe(3));
    [~,index] = min(abs(t - tframe));
    if t(index) >= tframe && index ~= 1, index = index - 1; end
    active = switches(index,:) & switches(index+1,:);
    set(splot, 'XData', vertices(active',1)*1.02, ...
               'YData', vertices(active',2)*1.02, ...
               'ZData', vertices(active',3)*1.02);
    drawnow;
    
    if record
        this_frame = getframe(fig);
        writeVideo(vidObj,this_frame);
    end
end
if record, close(vidObj); end

return
%%
dist = [];
for i=1:length(Eros.plates)
    r1 = Eros.vertices(Eros.plates(i,1)+1);
    r2 = Eros.vertices(Eros.plates(i,2)+1);
    r3 = Eros.vertices(Eros.plates(i,3)+1);
    dist(end+1) = norm(r1-r2);
    dist(end+1) = norm(r2-r3);
    dist(end+1) = norm(r3-r1);
end
mean_separation = mean(dist)*1e3
max_separation = max(dist)*1e3