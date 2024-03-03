function PaperAnimation
sim = load('sim_data_extended');

f = figure(1); clf;
scale = 1e3;
[~, vertices, indices] = make_die;
r = norm(mean(vertices(:,indices(:,1)),2));
rho = 1e4;
vertices = vertices*rho/r;
xc = zeros(size(sim.x));
for i=1:length(sim.x)
    xc(:,i) = get_center(sim.t(i));
end
vertices_pt = vertices + xc(1:3,1);
icosa = trimesh(indices',vertices_pt(1,:)/scale,vertices_pt(2,:)/scale,vertices_pt(3,:)/scale,'EdgeColor','k','FaceColor',[0.4;0.4;0.4],'FaceAlpha',0); hold on;

px = plot3(sim.x(1,1)/scale, sim.x(2,1)/scale, sim.x(3,1)/scale, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
pxhat = plot3(sim.xhat(1,1)/scale, sim.xhat(2,1)/scale, sim.xhat(3,1)/scale, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
[s1, s2, s3] = sphere(20);
fxhat = surf((s1*sim.rhohat(1,1) + sim.xhat(1,1))/scale,...
    (s2*sim.rhohat(1,1) + sim.xhat(2,1))/scale,...
    (s3*sim.rhohat(1,1) + sim.xhat(3,1))/scale, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeAlpha', 0);

axis equal;
xlabel 'x (km)'; ylabel 'y (km)'; zlabel 'z (km)';
[arrow1, arrow2] = draw_arrow(sim.x(1:3,1)/1e3, sim.u(:,1000)/norm(sim.u(:,1000)));
legend([px pxhat],{'Real','Estimate'});
set(f, 'Position', [50 100 1150 850]);
set(gca, 'FontSize', 13);
view(-200,6)
rho_min = 200;
t_info = title(['Time = ' num2str(0, '%.2f') ' hrs']);

moviename = 'Movie';
frame_rate = 40;
record = 1;
vidObj = VideoWriter(moviename);
vidObj.FrameRate = frame_rate;
if record, open(vidObj); end

for i=1:6:14401
    set(px, 'XData', sim.x(1,i)/scale, 'YData', sim.x(2,i)/scale, 'ZData', sim.x(3,i)/scale);
    vertices_pt = vertices + xc(1:3,i);
    set(icosa, 'Vertices', vertices_pt'/scale);
    set(fxhat, 'XData', (s1*max(rho_min, sim.rhohat(1,i)) + sim.xhat(1,i))/scale, ...
               'YData', (s2*max(rho_min, sim.rhohat(1,i)) + sim.xhat(2,i))/scale, ...
               'ZData', (s3*max(rho_min, sim.rhohat(1,i)) + sim.xhat(3,i))/scale);
    axis([-11 11 -11 11 -11 11] + [xc(1,i) xc(1,i) xc(2,i) xc(2,i) xc(3,i) xc(3,i)]/scale);
    set(t_info, 'String', ['Time = ' num2str(sim.t(i)/3600, '%.2f') ' hrs']);
    view(sim.t(i)/3600/7*20, 6);
    if norm(sim.u(:,i)) > 5e-5
        update_arrow(arrow1, arrow2, sim.x(1:3,i)/1e3, sim.u(:,i)/norm(sim.u(:,i)));
        set(arrow1, 'Visible', 'on');
        set(arrow2, 'Visible', 'on');
    else
        set(arrow1, 'Visible', 'off');
        set(arrow2, 'Visible', 'off');
    end
    drawnow;
    
    if record
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end

if record, close(vidObj); end

end

function [c, t] = draw_arrow(location, direction)
pts = get_arrow(location, direction);
x = reshape(pts(1,:), [2 21]);
y = reshape(pts(2,:), [2 21]);
z = reshape(pts(3,:), [2 21]);
c = surf(x, y, z, 'FaceColor', [0.6, 0.3 0.3], 'EdgeAlpha', 0);
x = reshape(pts(4,:), [2 21]);
y = reshape(pts(5,:), [2 21]);
z = reshape(pts(6,:), [2 21]);
t = surf(x, y, z, 'FaceColor', [0.6, 0.3 0.3], 'EdgeAlpha', 0);
end

function update_arrow(pc, pt, location, direction)
pts = get_arrow(location, direction);
x = reshape(pts(1,:), [2 21]);
y = reshape(pts(2,:), [2 21]);
z = reshape(pts(3,:), [2 21]);
set(pc, 'XData', x, 'YData', y, 'ZData', z);
x = reshape(pts(4,:), [2 21]);
y = reshape(pts(5,:), [2 21]);
z = reshape(pts(6,:), [2 21]);
set(pt, 'XData', x, 'YData', y, 'ZData', z);
end

function out = get_arrow(location, direction)
[t1, t2, t3] = cylinder(1,20);
[c1, c2, c3] = cylinder(1,20);
t1(2,:) = 0; t2(2,:) = 0;
tip_width = 0.4;
tip_height = 0.6;
base_height = 2;
base_width = 0.2;

t3(2,:) = 0;
t3(1,:) = -tip_height;
t2 = t2*tip_width;
t1 = t1*tip_width;

c3(2,:) = -tip_height;
c3(1,:) = -tip_height-base_height;
c2 = c2*base_width;
c1 = c1*base_width;

axis = cross([0;0;1], direction);
angle = acos(direction(3)/norm(direction));
q = [cos(angle/2); sin(angle/2)*axis/norm(axis)];

c_pts = [c1(:)'; c2(:)'; c3(:)'];
t_pts = [t1(:)'; t2(:)'; t3(:)'];
c_pts = RotQ(c_pts, q);
t_pts = RotQ(t_pts, q);

out = [c_pts + location;
       t_pts + location];
end
