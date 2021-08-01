load Results.mat
n = length(data);
load InData/Eros_Shape.mat;
t0 = 0;
t1 = 600;
video_duration = 60;
frame_rate = 50;
nframes = video_duration*frame_rate;

try
    load OutData/Results_Points.mat;
catch
    inputs = load('InData/Eros_Points.mat');
    points = [];
    for i=1:length(data)
        points{i} = [];
        for j=1:size(data(i).rc_all,1)
            site = data(i).rc_all(j,:);
            [~,index] = min(vecnorm(inputs.points - site,2,2));
            points{i}(j) = index;
        end
        waitbar(i/length(data));
    end
    save('OutData/Results_Points.mat', 'points');
end

h = 1-vecnorm([data.drc_alg]',2,2);
hdot = -dot([data.drc_alg], x(:,4:6)')'./(1-h);
H = h + max(hdot,0).^2/(u_max*2);
u = [data.u]'*1e3;
u_max = u_max*1e3;

%%

f1 = figure(1); clf;
set(gcf, 'Position', [100 100 1000 600]);
trisurf(plates+1, vertices(:,1), vertices(:,2), vertices(:,3), ...
    'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1, 'EdgeAlpha', 0.7);
hold on;
PlatesPlotHandle = trisurf([], Eros.vertices(:,1), Eros.vertices(:,2), Eros.vertices(:,3), ...
    'FaceColor', [1 0 0], 'FaceAlpha', 1, 'EdgeAlpha', 0);
PositionPlotHandle = plot3(x(1,1), x(1,2), x(1,3), 'b', 'LineWidth', 2);
CurrentPositionHandle = plot3(x(1,1), x(1,2), x(1,3), 'bo', 'MarkerFaceColor', [0;0;1]);
GroundPlotHandle = plot3(nan, nan, nan, 'go', 'MarkerFaceColor', [0;1;0]);
AxesHandle = gca;
xlabel 'x (km)'; ylabel 'y (km)'; zlabel 'z (km)'; hold on;
axis equal;
az = 176.1; el = 54.8;
view(az,el);

%%
f2 = figure(2); clf;
set(gcf, 'Position', [100 700 1000 250]);
subplot(1,2,1);
hPlotHandle = plot(nan, nan); hold on;
plot([0, 600], [0 0], 'r--');
xlabel 'Time (s)'; ylabel 'H (km)'; title 'ZCBF Condition';
bds{1} = gca;
subplot(1,2,2);
HPlotHandle = plot(nan, nan); hold on;
plot([0, 600], [0 0], 'r--');
xlabel 'Time (s)'; ylabel 'h (km)'; title 'Safety Condition';
bds{2} = gca;

tmin = 50;
f3 = figure(3); clf;
set(gcf, 'Position', [1100 100 500 850]);
subplot(3,1,1);
uxPlotHandle = plot(nan, nan); hold on;
plot([0 600], [-u_max, -u_max], 'r--');
plot([0 600], [u_max, u_max], 'r--');
ylabel 'u_x (m/s^2)'; title 'Control Input'
axis([0 tmin -u_max u_max]*1.2);
bds{3} = gca;

subplot(3,1,2);
uyPlotHandle = plot(nan, nan); hold on;
plot([0 600], [-u_max, -u_max], 'r--');
plot([0 600], [u_max, u_max], 'r--');
ylabel 'u_y (m/s^2)';
axis([0 tmin -u_max u_max]*1.2);
bds{4} = gca;

subplot(3,1,3);
uzPlotHandle = plot(nan, nan); hold on;
plot([0 600], [-u_max, -u_max], 'r--');
plot([0 600], [u_max, u_max], 'r--');
ylabel 'u_z (m/s^2)';
xlabel 'Time (s)';
axis([0 tmin -u_max u_max]*1.2);
bds{5} = gca;

moviename = 'Movie';
vidObj = VideoWriter(moviename);
vidObj.FrameRate = frame_rate;
open(vidObj);
this_frame = struct('cdata', [], 'colormap', []);

for frame=1:nframes
    tframe = frame/nframes*(t1 - t0) + t0;
    [~, index_max] = max(tframe < t);
    indices = 1:index_max;
    xframe = interp1(t, x, tframe);
    set(PositionPlotHandle, 'XData', x(indices,1), 'YData', x(indices,2), 'ZData', x(indices,3));
    set(CurrentPositionHandle, 'XData', xframe(1), 'YData', xframe(2), 'ZData', xframe(3));
    set(PlatesPlotHandle, 'Faces', Eros.plates(unique([points{index_max}, points{min(index_max+1,n)}]),:)+1);
    
    vec = xframe(1:3)/norm(xframe(1:3));
    az = atan2(vec(1),-vec(2))*180/pi+30;
    el = atan2(vec(3),sqrt(vec(1)^2+vec(2)^2))*180/pi+30;
    set(AxesHandle, 'View', [az el]);
    
    [~, ground_vec] = GroundMotion(xframe');
    set(GroundPlotHandle, 'XData', ground_vec(1), 'YData', ground_vec(2), 'ZData', ground_vec(3));
    drawnow;
    
    set(hPlotHandle, 'XData', t(indices), 'YData', h(indices));
    set(HPlotHandle, 'XData', t(indices), 'YData', H(indices));
    set(uxPlotHandle, 'XData', t(indices), 'YData', u(indices,1));
    set(uyPlotHandle, 'XData', t(indices), 'YData', u(indices,2));
    set(uzPlotHandle, 'XData', t(indices), 'YData', u(indices,3));
    
    for i=1:5
        set(bds{i}, 'XLim', [0 tmin*ceil(tframe/tmin)]);
    end
    
    
    frame_LeftTop = getframe(f1);
    frame_LeftBottom = getframe(f2);
    frame_Right = getframe(f3);
    this_frame.cdata = [[frame_LeftTop.cdata; frame_LeftBottom.cdata], frame_Right.cdata];
    writeVideo(vidObj,this_frame);
    
    waitbar(frame/nframes);
end

close(vidObj);