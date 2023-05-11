function Animation_Publish
sims = {};
sims{1} = load('Runs/linear30.mat');
sims{2} = load('Runs/linear45.mat');
sims{3} = load('Runs/linear60.mat');
sims{4} = load('Runs/nonlin30.mat');
sims{5} = load('Runs/nonlin45.mat');
sims{6} = load('Runs/nonlin60.mat');
sims{7} = load('Runs/nonlin180.mat');
sims{8} = load('Runs/nonlin300.mat');
sims{9} = load('Runs/nonlin420.mat');

f = figure(11); clf;
theta = linspace(0, 2*pi, 100);
cx = cos(theta);
cy = sin(theta);
rho = 0.85;
[xc, theta] = get_center(0);
center = xc(1:2)/1e3 + 6*[sin(theta); -cos(theta)];
obs = obstacle_location(1,0); o1 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r'); hold on;
obs = obstacle_location(2,0); o2 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(3,0); o3 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(4,0); o4 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
dock_width = 100;
o5 = fill(xc(1)/1e3 + dock_width*[-1, 1, 1, -1, -1], xc(2)/1e3 + dock_width*[0, 0, 1, 1, 0], 'r');
o_array = [o1; o2; o3; o4];
xlabel('x_1 (km)', 'FontSize', 14); ylabel('x_2 (km)', 'FontSize', 14);
axis equal;
axis([center(1)-8, center(1)+8, center(2)-8, center(2)+8]);
pc = plot(xc(1), xc(2), 'k.', 'MarkerSize', 20);
set(gcf, 'Position', [1700 300 760 600]);
lt = title([' ' newline], 'FontSize', 16);
lt = annotation('textbox',[0,0.894,0.9,0.1],'String','Impulsive Control Demonstration','FitBoxToText','on',...
    'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','LineStyle','none');
pt = annotation('textbox',[0,0.85,0.9,0.1],'String','t = 0.00 minutes','FitBoxToText','on',...
    'FontSize',13,'FontWeight','bold','HorizontalAlignment','center','LineStyle','none');

color1 = [0; 0.7; 0];
color2 = [0.2; 0.2; 1];
color3 = [0.6; 0; 0.6];
px = {};
px{1} = plot(sims{1}.x(1,1), sims{1}.x(2,1), 'ko', 'MarkerFaceColor', color1);
px{2} = plot(sims{2}.x(1,1), sims{2}.x(2,1), 'ko', 'MarkerFaceColor', color1+0.25);
px{3} = plot(sims{3}.x(1,1), sims{3}.x(2,1), 'ko', 'MarkerFaceColor', min(color1+0.5,1));
px{4} = plot(sims{4}.x(1,1), sims{4}.x(2,1), 'ko', 'MarkerFaceColor', color2);
px{5} = plot(sims{5}.x(1,1), sims{5}.x(2,1), 'ko', 'MarkerFaceColor', min(color2+0.25,1));
px{6} = plot(sims{6}.x(1,1), sims{6}.x(2,1), 'ko', 'MarkerFaceColor', min(color2+0.5,1));
px{7} = plot(sims{7}.x(1,1), sims{7}.x(2,1), 'ko', 'MarkerFaceColor', color3);
px{8} = plot(sims{8}.x(1,1), sims{8}.x(2,1), 'ko', 'MarkerFaceColor', min(color3+0.25,1));
px{9} = plot(sims{9}.x(1,1), sims{9}.x(2,1), 'ko', 'MarkerFaceColor', min(color3+0.5,1));

arrows = {};
arrows{1} = draw_arrow(sims{1}.x(1:2,1), sims{1}.x(3:4,1), '-', color1, 2);
arrows{2} = draw_arrow(sims{2}.x(1:2,1), sims{2}.x(3:4,1), '--', color1, 2);
arrows{3} = draw_arrow(sims{3}.x(1:2,1), sims{3}.x(3:4,1), ':', color1, 3);
arrows{4} = draw_arrow(sims{4}.x(1:2,1), sims{4}.x(3:4,1), '-', color2, 2);
arrows{5} = draw_arrow(sims{5}.x(1:2,1), sims{5}.x(3:4,1), '--', color2, 2);
arrows{6} = draw_arrow(sims{6}.x(1:2,1), sims{6}.x(3:4,1), ':', color2, 3);
arrows{7} = draw_arrow(sims{7}.x(1:2,1), sims{7}.x(3:4,1), '-', color3, 2);
arrows{8} = draw_arrow(sims{8}.x(1:2,1), sims{8}.x(3:4,1), '--', color3, 2);
arrows{9} = draw_arrow(sims{9}.x(1:2,1), sims{9}.x(3:4,1), ':', color3, 3);
arrow_times = -100*ones(1,9);
arrow_dwell_time = 15;

curr_indices = ones(1,9);
curr_u = ones(2,9);
last_u = ones(2,9);

run_fast = 0;
if ~run_fast % I am not sure why the presence of a legend slows down the animator so much
le = legend([px{:}], {'$$\psi_h$$, 30 s', '$$\psi_h$$, 45 s', '$$\psi_h$$, 60 s', ...
            '$$\psi_h^*$$, 30 s', '$$\psi_h^*$$, 45 s', '$$\psi_h^*$$, 60 s', ...
            '$$\psi_h^*$$, 180 s', '$$\psi_h^*$$, 300 s', '$$\psi_h^*$$, 420 s'}, ...
       'Location', 'NorthEast', 'FontSize', 13, 'interpreter', 'latex');
ax = gca;
set(ax, 'Position', ax.Position - [0.07, 0, 0, 0]);
set(le, 'Position', [0.79   0.578   0.157046917694350   0.325333326896032])
end

moviename = 'Movie';
frame_rate = 40;
record = 1;
vidObj = VideoWriter(moviename);
vidObj.FrameRate = frame_rate;
if record, open(vidObj); end

speed = 1;
t = 0;
for i=1:5000
    for j=1:9
        if sims{j}.t(curr_indices(j)+1) == sims{j}.t(curr_indices(j))
            curr_u(:,j) = sims{j}.u(:,curr_indices(j));
            last_u(:,j) = curr_u(:,j);
            curr_indices(j) = curr_indices(j)+1;
        else
            curr_u(:,j) = [0;0];
        end
        if sims{j}.t(curr_indices(j)) ~= t && ~isnan(sims{j}.t(curr_indices(j)))
            disp('Timing Mistmatch');
        end
    end
    
    for j=1:4
        obs = obstacle_location(j,t);
        set(o_array(j), 'Vertices', [obs(1)/1e3+cx(:)*rho, obs(2)/1e3+cy(:)*rho]);
    end
    [xc, theta] = get_center(t);
    rectangle = [20*[-1, 1, 1, -1, -1]; 0, 0.15, 10, 10, 0.15];
    rectangle = [cos(theta), -sin(theta); sin(theta), cos(theta)]*rectangle;
    set(o5, 'Vertices', [xc(1), xc(2)]/1e3 + rectangle');
    set(pc, 'XData', xc(1)/1e3, 'YData', xc(2)/1e3);
    center = xc(1:2)/1e3 + 6*[sin(theta); -cos(theta)];
    axis([center(1)-8, center(1)+8, center(2)-8, center(2)+8]);
    set(pt, 'String', ['t = ' num2str(t/60, '%.2f') ' minutes']);
    
    center_text = round(center(1)/2)*2;
    set(ax, 'XTick', center_text+(-6:2:6));
    center_text = round(center(2)/2)*2;
    set(ax, 'YTick', center_text+(-6:2:6));
    
    for j=1:9
        k = curr_indices(j);
        if k<= length(sims{j}.t) && ~isnan(sims{j}.t(k))
            set(px{j}, 'XData', sims{j}.x(1,k)/1e3, 'YData', sims{j}.x(2,k)/1e3);
            pts = get_arrow(sims{j}.x(1:2,k)/1e3, last_u(:,j));
            set(arrows{j}, 'XData', pts(1,:), 'YData', pts(2,:));
            if norm(curr_u(:,j)) ~= 0
                set(arrows{j}, 'Visible', 'on');
                arrow_times(j) = t;
            end
        else
            disp(['Reached the end of simulation ' num2str(j)]);
        end
        if t > arrow_times(j) + arrow_dwell_time*speed
            set(arrows{j}, 'Visible', 'off');
        end
    end
    
    if t > 750
        speed = 2;
    end
    if t > 2250
        speed = 4;
    end
    if mod(t,speed)==0
        drawnow;
        if record
            frame = getframe(f);
            writeVideo(vidObj,frame);
        end
    end
    t = t+1;
    curr_indices = curr_indices+1;
end

if record, close(vidObj); end
end

function p = draw_arrow(location, direction, style, color, width)
pts = get_arrow(location, direction);
p = plot(pts(1,:), pts(2,:), style, 'LineWidth', width, 'Color', color);
end

function pts = get_arrow(location, direction)
tip_length = 0.3;
base_length = 1.5;
x = [-tip_length, 0, tip_length, NaN, 0, 0];
y = [-1.5*tip_length, 0, -1.5*tip_length, NaN, 0, -base_length];
pts = [x; y];
theta = atan2(direction(2), direction(1));
pts = [cos(theta), -sin(theta); sin(theta), cos(theta)]*pts + [location(1); location(2)];
end