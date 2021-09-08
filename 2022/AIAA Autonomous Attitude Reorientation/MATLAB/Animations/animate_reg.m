%% Initialization
fig1 = figure(1); clf;
cubesat = loadmesh; hold on;

sphere_density = 300;

constants = load('../parameters.mat');
z = [0;0;-1];
y = [0;1;0];
theta = acos(dot(constants.p1, z));
v = cross(constants.p1, z);
dq1 = [cos(theta/2); v/norm(v)*sin(theta/2)]; % active rotation form p1 to z

p2 = RotQ(constants.p2, dq1);
p2 = p2 - dot(p2,z)*z; p2 = p2/norm(p2);
theta = acos(dot(p2, y));
v = cross(p2, y);
dq2 = [cos(theta/2); v/norm(v)*sin(theta/2)]; % active rotation from p2 projection to z

dq = QxQ(dq2, dq1); % rotation to better align p1 and p2 with the cubesat model
    % p1 and p2 were chosen arbitrarily before the cubesat model was actually made, but
    % for this visualization, we want them to be a little more meaningful.

[xs,ys,zs] = sphere(40);
sphere_radius = 55;
surf(xs*sphere_radius,ys*sphere_radius,zs*sphere_radius, 'FaceAlpha', 0, 'EdgeAlpha', 0.2, 'EdgeColor', [0.4;0.4;0.4]);
[xs,ys,zs] = sphere(sphere_density); 

vector = RotQ(constants.p1, dq); % minus z axis, as we desire
[x,y,z] = cylinder(1,40);
height = 12;
angle = acos(constants.ctheta);
radius = height*tan(angle);
offset = 15;
instrument1 = surf(x*radius.*[0;1], y*radius.*[0;1]-6, -z*height-offset, ...
    'FaceColor', [0; 0.6; 0], 'EdgeAlpha', 0, 'FaceAlpha', 0.7);
indices = (xs*vector(1) + ys*vector(2) + zs*vector(3)) >= cos(angle);
x1 = xs*sphere_radius; x1(~indices) = nan;
y1 = ys*sphere_radius; y1(~indices) = nan;
z1 = zs*sphere_radius; z1(~indices) = nan;
instrument_sphere1 = surf(x1, y1, z1, 'FaceColor', [0; 0.6; 0], 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
xc = -1:0.02:1;
yc = sqrt(1-xc.^2);
xc = [xc, fliplr(xc)]; yc = [yc, -yc];
x1 = xc*sin(angle)*sphere_radius;
y1 = yc*sin(angle)*sphere_radius;
z1 = -ones(size(xc))*cos(angle)*sphere_radius;
instrument_circle1 = plot3(x1, y1, z1, 'Color', [0; 0.6; 0], 'LineWidth', 4);

vector = RotQ(constants.p2, dq); % approximately -y axis
height = 12;
angle = acos(constants.ctheta);
radius = height*tan(angle);
x1 = x*radius.*[0;1];
y1 = y*radius.*[0;1];
z1 = z*height;
[x2,y2,z2] = rotate_coordinates(get_rotation([0;0;1],vector), x1, y1, z1);
instrument2 = surf(x2, y2+10, z2+11, ...
    'FaceColor', [0.2; 0.2; 1], 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
indices = (xs*vector(1) + ys*vector(2) + zs*vector(3)) >= cos(angle);
x1 = xs*sphere_radius; x1(~indices) = nan;
y1 = ys*sphere_radius; y1(~indices) = nan;
z1 = zs*sphere_radius; z1(~indices) = nan;
instrument_sphere2 = surf(x1, y1, z1, 'FaceColor', [0.2; 0.2; 1], 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
x1 = xc*sin(angle)*sphere_radius;
y1 = yc*sin(angle)*sphere_radius;
z1 = -ones(size(xc))*cos(angle)*sphere_radius;
[x2, y2, z2] = rotate_coordinates(get_rotation([0;0;-1],vector), x1, y1, z1);
instrument_circle2 = plot3(x2, y2, z2, 'Color', [0.2; 0.2; 1], 'LineWidth', 4);

vector = [1;0;0]; % get_s(0);
top = 55;
sun_vector = plot3([0, top], [0, 0], [0, 0], 'Color', [.9; 1; 0], 'LineWidth', 3);
radius = 1;
height = 2.5;
x1 = x*radius.*[0;1];
y1 = y*radius.*[0;1];
z1 = z*height;
[x2,y2,z2] = rotate_coordinates(get_rotation([0;0;1],vector), x1, y1, z1);
sunarrow = surf(top+1-x2, y2, z2, ...
    'FaceColor', [0.9; 1; 0], 'EdgeAlpha', 0, 'FaceAlpha', 1);
sun_sphere = [];%plot3(sphere_radius, 0, 0, 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 6);

s_target = [0; -0.707106781186547; -0.707106781186547];
R = get_rotation([0;0;1], s_target);
top = 55;
target_vector = plot3([0, s_target(1)*top], [0, s_target(2)*top], [0, s_target(3)*top],...
    'Color', [0; 0.9; 0], 'LineWidth', 3);
radius = 1;
height = 2.5;
[cx, cy, cz] = rotate_coordinates(R, x*radius.*[0; 1], y*radius.*[0;1], top+1-z*height);
target_arrow = surf(cx, cy, cz, 'FaceColor', [0; 0.9; 0], 'EdgeAlpha', 0, 'FaceAlpha', 1);
target_sphere = [];%plot3(sphere_radius*s_target(1), sphere_radius*s_target(2), sphere_radius*s_target(3),...
%     'o', 'Color', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', 6);

spacecraft1 = struct('meshes', cubesat, ...
    'surfaces', [instrument1; instrument2; instrument_sphere1; instrument_sphere2], ...
    'lines', [instrument_circle1; instrument_circle2], 'fixed_lines', [], ...
    'fixed_surfaces', [], ...
    'fixed_position', [0; 0; 0], 'R', eye(3));
    % This code is adapted from earlier code, so some features are no actually used.
spacecraft1 = rotate_spacecraft(QtoR(QConj(dq)), spacecraft1);

set(fig1, 'Position', [100 100 800 700]);
axis off;
a = gca;
set(a, 'XLim', [-60 60])
set(a, 'YLim', [-60 60])
set(a, 'ZLim', [-60 70])
axis equal;
set(fig1, 'color', 'k');
view(52,2);


%% Our approach
fig2 = figure(2); clf;
cubesat = loadmesh; hold on;

[xs,ys,zs] = sphere(40);
surf(xs*sphere_radius,ys*sphere_radius,zs*sphere_radius, 'FaceAlpha', 0, 'EdgeAlpha', 0.2, 'EdgeColor', [0.4;0.4;0.4]);
[xs,ys,zs] = sphere(sphere_density); 

vector = RotQ(constants.p1, dq); % minus z axis, as we desire
[x,y,z] = cylinder(1,40);
height = 12;
angle = acos(constants.ctheta);
radius = height*tan(angle);
offset = 15;
instrument1 = surf(x*radius.*[0;1], y*radius.*[0;1]-6, -z*height-offset, ...
    'FaceColor', [0; 0.6; 0], 'EdgeAlpha', 0, 'FaceAlpha', 0.7);
indices = (xs*vector(1) + ys*vector(2) + zs*vector(3)) >= cos(angle);
x1 = xs*sphere_radius; x1(~indices) = nan;
y1 = ys*sphere_radius; y1(~indices) = nan;
z1 = zs*sphere_radius; z1(~indices) = nan;
instrument_sphere1 = surf(x1, y1, z1, 'FaceColor', [0; 0.6; 0], 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
xc = -1:0.02:1;
yc = sqrt(1-xc.^2);
xc = [xc, fliplr(xc)]; yc = [yc, -yc];
x1 = xc*sin(angle)*sphere_radius;
y1 = yc*sin(angle)*sphere_radius;
z1 = -ones(size(xc))*cos(angle)*sphere_radius;
instrument_circle1 = plot3(x1, y1, z1, 'Color', [0; 0.6; 0], 'LineWidth', 4);

vector = RotQ(constants.p2, dq); % approximately -y axis
height = 12;
angle = acos(constants.ctheta);
radius = height*tan(angle);
x1 = x*radius.*[0;1];
y1 = y*radius.*[0;1];
z1 = z*height;
[x2,y2,z2] = rotate_coordinates(get_rotation([0;0;1],vector), x1, y1, z1);
instrument2 = surf(x2, y2+10, z2+11, ...
    'FaceColor', [0.2; 0.2; 1], 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
indices = (xs*vector(1) + ys*vector(2) + zs*vector(3)) >= cos(angle);
x1 = xs*sphere_radius; x1(~indices) = nan;
y1 = ys*sphere_radius; y1(~indices) = nan;
z1 = zs*sphere_radius; z1(~indices) = nan;
instrument_sphere2 = surf(x1, y1, z1, 'FaceColor', [0.2; 0.2; 1], 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
x1 = xc*sin(angle)*sphere_radius;
y1 = yc*sin(angle)*sphere_radius;
z1 = -ones(size(xc))*cos(angle)*sphere_radius;
[x2, y2, z2] = rotate_coordinates(get_rotation([0;0;-1],vector), x1, y1, z1);
instrument_circle2 = plot3(x2, y2, z2, 'Color', [0.2; 0.2; 1], 'LineWidth', 4);

vector = [1;0;0]; % get_s(0);
top = 55;
sun_vector = plot3([0, top], [0, 0], [0, 0], 'Color', [.9; 1; 0], 'LineWidth', 3);
radius = 1;
height = 2.5;
x1 = x*radius.*[0;1];
y1 = y*radius.*[0;1];
z1 = z*height;
[x2,y2,z2] = rotate_coordinates(get_rotation([0;0;1],vector), x1, y1, z1);
sunarrow = surf(top+1-x2, y2, z2, ...
    'FaceColor', [0.9; 1; 0], 'EdgeAlpha', 0, 'FaceAlpha', 1);
sun_sphere = [];%plot3(sphere_radius, 0, 0, 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 6);

s_target = [0; -0.707106781186547; -0.707106781186547];
R = get_rotation([0;0;1], s_target);
top = 55;
target_vector = plot3([0, s_target(1)*top], [0, s_target(2)*top], [0, s_target(3)*top],...
    'Color', [0; 0.9; 0], 'LineWidth', 3);
radius = 1;
height = 2.5;
[cx, cy, cz] = rotate_coordinates(R, x*radius.*[0; 1], y*radius.*[0;1], top+1-z*height);
target_arrow = surf(cx, cy, cz, 'FaceColor', [0; 0.9; 0], 'EdgeAlpha', 0, 'FaceAlpha', 1);
target_sphere = [];%plot3(sphere_radius*s_target(1), sphere_radius*s_target(2), sphere_radius*s_target(3),...
%     'o', 'Color', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', 6);

spacecraft2 = struct('meshes', cubesat, ...
    'surfaces', [instrument1; instrument2; instrument_sphere1; instrument_sphere2], ...
    'lines', [instrument_circle1; instrument_circle2], 'fixed_lines', [], ...
    'fixed_surfaces', [], ...
    'fixed_position', [0; 0; 0], 'R', eye(3));
    % This code is adapted from earlier code, so some features are no actually used.
spacecraft2 = rotate_spacecraft(QtoR(QConj(dq)), spacecraft2);

set(fig2, 'Position', [900 100 800 700]);
axis off;
a = gca;
set(a, 'XLim', [-60 60])
set(a, 'YLim', [-60 60])
set(a, 'ZLim', [-60 70])
axis equal;
set(fig2, 'color', 'k');
view(52,2);


%% Animation
results_left = load('../Results/Comparison.mat');
results_right = load('../Results/Nominal.mat');
t0 = 0;
t1 = results_left.SimDur;

video_duration = 60;
frame_rate = 30; %10; % for faster previews.
nframes = video_duration*frame_rate;

record = 0;

moviename = 'MovieNom2';
vidObj = VideoWriter(moviename);
vidObj.FrameRate = frame_rate;
if record, open(vidObj); end

%% Pre-Info
fig3 = figure(3); clf
set(fig3, 'Position', [100 100 1600 700]);
axis off;
a = gca;
set(a, 'XLim', [0 300])
set(a, 'YLim', [-60 60])
set(fig3, 'color', 'k');
text(-40, 36, ['"Autonomous Spacecraft Attitude Reorientation' newline ...
    '                Using Robust Sampled-Data' newline ...
    '                  Control Barrier Functions"'],...
    'Color', 'w', 'Fontsize', 48, 'FontWeight', 'Bold');
text(10, -14, 'Comparison of a Nominal Reorientation', 'Color', 'w', 'Fontsize', 44);
text(12.5, -15.5, '_______________________________', 'color', 'w', 'fontsize', 44);
text(39, -44, 'Joseph Breeden and Dimitra Panagou', 'Color', 'w', 'Fontsize', 35);

duration = 5;
fade = 1;
if record
    for i=1:duration*frame_rate
        this_frame = getframe(fig3);
        writeVideo(vidObj,this_frame);
    end
    fade_frames = round(fade*frame_rate);
    for i=1:fade_frames
        this_frame = getframe(fig3);
        new_im = this_frame.cdata;
        m = max(new_im(:));
        new_im = new_im - (m*i/fade_frames);
        this_frame.cdata = new_im;
        writeVideo(vidObj,this_frame);
    end
end

surf(a.XLim, a.YLim, ones(2), 'FaceColor', 'k');
view([0 0 1]);

text(-30,40,'The yellow vector is the current Sun vector.','Color',[1;1;0],'FontSize',40);
text(-30,28,'The sensitive instruments are the           and','color','w','fontsize',40); 
text(163.3,28,'green','color',[0;0.6;0],'fontsize',40)
text(223,28,'blue','color',[0.2;0.2;1],'fontsize',40)
text(250,28,'regions.','color','w','fontsize',40);

text(-30, 0, 'The spacecraft must reorient so the green instrument','color','w','fontsize',40);
text(-25, -12, 'points at the ','color','w','fontsize',40); 
text(49, -12, 'green arrow','color','g','fontsize',40); 
text(117, -12, ', while both instruments','color','w','fontsize',40);
text(-25, -24, 'avoid the Sun vector.','color','w','fontsize',40);

duration = 5;
fade = 1;
if record
    for i=1:duration*frame_rate
        this_frame = getframe(fig3);
        writeVideo(vidObj,this_frame);
    end
    fade_frames = round(fade*frame_rate);
    for i=1:fade_frames
        this_frame = getframe(fig3);
        new_im = this_frame.cdata;
        m = max(new_im(:));
        new_im = new_im - (m*i/fade_frames);
        this_frame.cdata = new_im;
        writeVideo(vidObj,this_frame);
    end
end

%% Spacecraft Mission
% Bring to front
figure(fig1); 
text(-54, 0, 68, 'Comparison', 'color', 'w', 'fontsize', 40);

figure(fig2);
text(-36, 0, 68, 'ZohCBF', 'color', 'w', 'fontsize', 40);

q0_left = results_left.q0;
spacecraft1 = rotate_spacecraft(QtoR(q0_left), spacecraft1);
q0_right = results_right.q0;
spacecraft2 = rotate_spacecraft(QtoR(q0_right), spacecraft2);

event1 = 0;
event2 = 0;
event3 = 0;
event4 = 0;
event5 = 0;
event6 = 0;
event7 = 0;
event8 = 0;
event9 = 0;
event10 = 0;
% nframes = 0; % this line is to let us skip the animation.
for frame=1:nframes
    t = t0 + (t1-t0)*frame/nframes;
    q1 = interp1(results_left.t, results_left.x(:,1:4), t); q1 = q1(:);
    R = QtoR(QxQ(q1, QConj(q0_left)));
    q0_left = q1;
    spacecraft1 = rotate_spacecraft(R, spacecraft1);
    
    q1 = interp1(results_right.t, results_right.x(:,1:4), t); q1 = q1(:);
    R = QtoR(QxQ(q1, QConj(q0_right)));
    q0_right = q1;
    spacecraft2 = rotate_spacecraft(R, spacecraft2);
    
    if t >= 20 && ~event1
        figure(fig1);
        h1 = text(-130, 0, -65, 'The comparison controller reoriented the','color','w','fontsize',24);
        h2 = text(-130, 0, -73, 'instrument above the sun vector.', 'color','w', 'fontsize',24);
        
        figure(fig2);
        h3 = text(-130, 0, -65, 'The ZohCBF controller reoriented the','color','w','fontsize',24);
        h4 = text(-130, 0, -73, 'instrument below the Sun vector.', 'color','w', 'fontsize',24);
        event1 = 1;
    end
    
    if t >= 56 && ~event2 % one minute of simtime is about ten seconds display time
        delete(h1);
        delete(h2);
        delete(h3);
        delete(h4);
        event2 = 1;
    end
    
    if t >= 60 && ~event3
        figure(fig2);
        h3 = text(-140, 0, -65, 'The ZohCBF controller also approched much','color','w','fontsize',24);
        h4 = text(-140, 0, -73, 'closer to the edge of the safe set.', 'color','w', 'fontsize',24);
        event3 = 1;
    end
    
    if t>=100 && ~event4
        delete(h3);
        delete(h4);
        event4 = 1;
    end
    
    if t>= 115 && ~event5
        figure(fig2);
        h3 = text(-140, 0, -65, 'Because of the chosen route, the ZohCBF','color','w','fontsize',24);
        h4 = text(-140, 0, -73, 'controller now must guide the second','color','w','fontsize',24);
        h5 = text(-140, 0, -81, 'instrument (blue) around the Sun vector.', 'color','w', 'fontsize',24);
        event5 = 1;
    end
    
    if t>=175 && ~event6
        delete(h3);
        delete(h4);
        delete(h5);
        event6 = 1;
    end
        
    if t>= 200 && ~event7
        figure(fig1);
        h1 = text(-140, 0, -65, 'Due to the tuning to meet the angular velocity','color','w','fontsize',24);
        h2 = text(-140, 0, -73, 'requirement, the comparison controller','color','w','fontsize',24);
        h3 = text(-140, 0, -81, 'converges more slowly.', 'color','w', 'fontsize',24);
        event7 = 1;
    end
    
    if t>= 250 && ~event8
        delete(h1);
        delete(h2);
        delete(h3);
        event8 = 1;
    end
    
    if t>= 300 && ~event9
        figure(fig1);
        h1 = text(-140, 0, -65, 'Both controllers accomplished the mission.','color','w','fontsize',24);
        figure(fig2);
        h3 = text(-140, 0, -65, 'Both controllers accomplished the mission.','color','w','fontsize',24);
        event9 = 1;
    end
        
    
    drawnow;
    
	if record
        left_frame = getframe(fig1);
        right_frame = getframe(fig2);
        this_frame = left_frame;
        this_frame.cdata = [left_frame.cdata, right_frame.cdata];
        writeVideo(vidObj,this_frame);
    end
    
    waitbar(frame/nframes);
end

if record, close(vidObj); end
