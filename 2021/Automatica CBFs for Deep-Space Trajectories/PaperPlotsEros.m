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
%%
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