function PaperAnimation
sim1 = load('data/sim_data120');
sim2 = load('data/sim_data240');
sim3 = load('data/sim_data360');

f = figure(11); clf;
theta = linspace(0, 2*pi, 100);
cx = cos(theta);
cy = sin(theta);
rho = 1;
xc = get_center(0);
obs = obstacle_location(1,0); o1 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r'); hold on;
obs = obstacle_location(2,0); o2 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(3,0); o3 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(4,0); o4 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(5,0); o5 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(6,0); o6 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
dock_width = 100;
o7 = fill(xc(1)/1e3 + dock_width*[-1, 1, 1, -1, -1], xc(2)/1e3 + dock_width*[0, 0, 1, 1, 0], 'r');
o_array = [o1; o2; o3; o4; o5; o6];
xlabel 'x_1 (km)'; ylabel 'x_2 (km)';
px1 = plot(sim1.x(1,1)/1e3, sim1.x(2,1)/1e3, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
px2 = plot(sim2.x(1,1)/1e3, sim2.x(2,1)/1e3, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
px3 = plot(sim3.x(1,1)/1e3, sim3.x(2,1)/1e3, 'ko', 'MarkerFaceColor', 'm', 'MarkerSize', 6);
pc = plot(xc(1)/1e3, xc(2)/1e3, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

arrow1 = draw_arrow(sim1.x(1:2,1)/1e3, sim1.x(3:4,1)/1e3, '--', 'b', 2);
arrow2 = draw_arrow(sim2.x(1:2,1)/1e3, sim2.x(3:4,1)/1e3, '--', 'g', 2);
arrow3 = draw_arrow(sim3.x(1:2,1)/1e3, sim3.x(3:4,1)/1e3, '--', 'm', 2);
arrow_time1 = -100;
arrow_time2 = -100;
arrow_time3 = -100;
arrow_dwell_time = 80;
last_u1 = [0;0];
last_u2 = [0;0];
last_u3 = [0;0];

rx1 = sim1.x(1,1) + cx*sim1.rhohat(1,1);
ry1 = sim1.x(2,1) + cy*sim1.rhohat(2,1);
f1 = fill(rx1, ry1, 'b', 'FaceAlpha', 1);
rx2 = sim2.x(1,1) + cx*sim2.rhohat(1,1);
ry2 = sim2.x(2,1) + cy*sim2.rhohat(2,1);
f2 = fill(rx2, ry2, 'g', 'FaceAlpha', 1);
rx3 = sim3.x(1,1) + cx*sim3.rhohat(1,1);
ry3 = sim3.x(2,1) + cy*sim3.rhohat(2,1);
f3 = fill(rx3, ry3, 'm', 'FaceAlpha', 1);

axis equal;
axis_size = [-20  20 -32 5];
axis([xc(1) xc(1) xc(2) xc(2)]/1e3 + axis_size);
axis_pts = [axis_size(1), axis_size(1), axis_size(2), axis_size(2); axis_size(3), axis_size(4), axis_size(3), axis_size(4)];
legend([px1 px2 px3], {'120 s', '240 s', '360 s'}, 'FontSize', 13, 'Location', 'SouthEast');

i1 = 1;
i2 = 1;
i3 = 1;

set(f, 'Position', [260 100 900 800])

moviename = 'Movie2';
frame_rate = 36;
record = 1;
vidObj = VideoWriter(moviename);
vidObj.FrameRate = frame_rate;
if record, open(vidObj); end

t_info = title('Time = 0 s');

%%
for t = [0:2:2500, 2504:6:5000]
    if t==2500
        arrow_dwell_time = arrow_dwell_time*2;
    end
    
    set(t_info, 'String', ['Time = ' num2str(t) ' s']);
    
    while sim1.t(i1) < t, i1 = i1 + 1; end
    while sim2.t(i2) < t, i2 = i2 + 1; end
    while sim3.t(i3) < t, i3 = i3 + 1; end
    
    xi1 = sim1.x(:,i1);
    xi2 = sim2.x(:,i2);
    xi3 = sim3.x(:,i3);
    for j=1:6
        obs = obstacle_location(j,t);
        set(o_array(j), 'Vertices', [obs(1)/1e3+cx(:)*rho, obs(2)/1e3+cy(:)*rho]);
    end
    [xc, ~] = get_center(t);
    rectangle = [-1, 1, 1, -1, -1; 0, 0, 1, 1, 0];
    theta = atan2(xc(4), xc(3)) - pi/2;
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    rectangle = R*rectangle;
    set(o7, 'Vertices', [xc(1), xc(2)]/1e3 + dock_width*rectangle');
%     set(px1, 'XData', xi1(1)/1e3, 'YData', xi1(2)/1e3);
%     set(px2, 'XData', xi2(1)/1e3, 'YData', xi2(2)/1e3);
%     set(px3, 'XData', xi3(1)/1e3, 'YData', xi3(2)/1e3);
    set(pc, 'XData', xc(1)/1e3, 'YData', xc(2)/1e3);
    curr_axis_pts = R*axis_pts;
    
    arrow_pts1 = get_arrow(sim1.x(1:2,i1)/1e3, last_u1);
    set(arrow1, 'XData', arrow_pts1(1,:), 'YData', arrow_pts1(2,:));
    if norm(sim1.u(:,i1)) ~= 0 && ~isnan(sim1.u(1,i1)), set(arrow1, 'Visible', 'on'); arrow_time1 = t; last_u1 = sim1.u(:,i1); end
    if t > arrow_time1 + arrow_dwell_time, set(arrow1, 'Visible', 'off'); end
    
    arrow_pts2 = get_arrow(sim2.x(1:2,i2)/1e3, last_u2);
    set(arrow2, 'XData', arrow_pts2(1,:), 'YData', arrow_pts2(2,:));
    if norm(sim2.u(:,i2)) ~= 0 && ~isnan(sim2.u(1,i2)), set(arrow2, 'Visible', 'on'); arrow_time2 = t; last_u2 = sim2.u(:,i2); end
    if t > arrow_time2 + arrow_dwell_time, set(arrow2, 'Visible', 'off'); end
    
    arrow_pts3 = get_arrow(sim3.x(1:2,i3)/1e3, last_u3);
    set(arrow3, 'XData', arrow_pts3(1,:), 'YData', arrow_pts3(2,:));
    if norm(sim3.u(:,i3)) ~= 0 && ~isnan(sim3.u(1,i3)), set(arrow3, 'Visible', 'on'); arrow_time3 = t; last_u3 = sim3.u(:,i3); end
    if t > arrow_time3 + arrow_dwell_time, set(arrow3, 'Visible', 'off'); end
    
    rx1 = sim1.xhat(1,i1) + cx*max(200, sim1.rhohat(1,i1));
    ry1 = sim1.xhat(2,i1) + cy*max(200, sim1.rhohat(1,i1));
    set(f1, 'XData', rx1/1e3, 'YData', ry1/1e3);
    rx2 = sim2.xhat(1,i2) + cx*max(200, sim2.rhohat(1,i2));
    ry2 = sim2.xhat(2,i2) + cy*max(200, sim2.rhohat(1,i2));
    set(f2, 'XData', rx2/1e3, 'YData', ry2/1e3);
    rx3 = sim3.xhat(1,i3) + cx*max(200, sim3.rhohat(1,i3));
    ry3 = sim3.xhat(2,i3) + cy*max(200, sim3.rhohat(1,i3));
    set(f3, 'XData', rx3/1e3, 'YData', ry3/1e3);
    
    axis_size = [-15  15 max(min(-15, min([sim1.x_lin(2,i1), sim2.x_lin(2,i2), sim3.x_lin(2,i3)]/1e3-2)), -32) 3];
    axis_pts = [axis_size(1), axis_size(1), axis_size(2), axis_size(2); axis_size(3), axis_size(4), axis_size(3), axis_size(4)];
    axis_size = [min(curr_axis_pts(1,:)), max(curr_axis_pts(1,:)), min(curr_axis_pts(2,:)), max(curr_axis_pts(2,:))];
    axis([xc(1), xc(1), xc(2), xc(2)]/1e3 + axis_size);
    drawnow;
%     pause(0.01);
    
    if record
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end

if record, close(vidObj); end

end

function p = draw_arrow(location, direction, style, color, width)
pts = get_arrow(location, direction);
p = plot(pts(1,:), pts(2,:), style, 'LineWidth', width, 'Color', color);
end

function pts = get_arrow(location, direction)
tip_length = 0.4*2;
base_length = 1.8*3;
x = [-tip_length, 0, tip_length, NaN, 0, 0];
y = [-1.5*tip_length, 0, -1.5*tip_length, NaN, 0, -base_length];
pts = [x; y];
theta = atan2(direction(2), direction(1));
pts = [cos(theta), -sin(theta); sin(theta), cos(theta)]*pts + [location(1); location(2)];
end