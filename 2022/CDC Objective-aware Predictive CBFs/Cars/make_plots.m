function make_plots
data1 = load('data_OPCBF_parallel');
data2 = load('data_ECBF_parallel');
data3 = load('data_MPC_parallel');

data4 = load('data_OPCBF_left');
data5 = load('data_ECBF_left');
data6 = load('data_MPC_left');

%% Control Input Plots
f = figure(2); clf;
plot(data1.t(1:end-1), data1.u(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(data2.t(1:end-1), data2.u(1,:), 'g', 'LineWidth', 1.5);
plot(data3.t(1:end-1), data3.u(1,:), 'r', 'LineWidth', 1.5);
plot(data1.t(1:end-1), data1.u(2,:), 'b--', 'LineWidth', 1.5);
plot(data2.t(1:end-1), data2.u(2,:), 'g--', 'LineWidth', 1.5);
plot(data3.t(1:end-1), data3.u(2,:), 'r--', 'LineWidth', 1.5);
xlabel 'Time (s)';
ylabel 'u (m/s^2)';
set(f, 'Position', [2200 1100 560 180]);
axis([0 8 -5 2.5]);
legend({'OPCBF', 'ECBF', 'MPC'},'Location','SouthEast','FontSize',12);

f = figure(3); clf;
plot(data4.t(1:end-1), data4.u(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(data5.t(1:end-1), data5.u(1,:), 'g', 'LineWidth', 1.5);
plot(data6.t(1:end-1), data6.u(1,:), 'r', 'LineWidth', 1.5);
plot(data4.t(1:end-1), data4.u(2,:), 'b--', 'LineWidth', 1.5);
plot(data5.t(1:end-1), data5.u(2,:), 'g--', 'LineWidth', 1.5);
plot(data6.t(1:end-1), data6.u(2,:), 'r--', 'LineWidth', 1.5);
xlabel 'Time (s)';
ylabel 'u (m/s^2)';
set(f, 'Position', [2200 800 560 180]);
axis([0 8 -5 2.5]);
legend({'OPCBF', 'ECBF', 'MPC'},'Location','SouthEast','FontSize',12);

%% Control Input Plots
f = figure(2); clf;
plot(data1.t(1:end-1), data1.u(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(data2.t(1:end-1), data2.u(1,:), 'g', 'LineWidth', 1.5);
plot(data3.t(1:end-1), data3.u(1,:), 'r', 'LineWidth', 1.5);
plot(data1.t(1:end-1), data1.u(2,:), 'b--', 'LineWidth', 1.5);
plot(data2.t(1:end-1), data2.u(2,:), 'g--', 'LineWidth', 1.5);
plot(data3.t(1:end-1), data3.u(2,:), 'r--', 'LineWidth', 1.5);
xlabel 'Time (s)';
ylabel 'u (m/s^2)';
set(f, 'Position', [2200 1100 560 180]);
axis([0 8 -5 2.5]);
legend({'OPCBF', 'ECBF', 'MPC'},'Location','SouthEast','FontSize',12);

f = figure(3); clf;
plot(data4.t(1:end-1), data4.u(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(data5.t(1:end-1), data5.u(1,:), 'g', 'LineWidth', 1.5);
plot(data6.t(1:end-1), data6.u(1,:), 'r', 'LineWidth', 1.5);
plot(data4.t(1:end-1), data4.u(2,:), 'b--', 'LineWidth', 1.5);
plot(data5.t(1:end-1), data5.u(2,:), 'g--', 'LineWidth', 1.5);
plot(data6.t(1:end-1), data6.u(2,:), 'r--', 'LineWidth', 1.5);
xlabel 'Time (s)';
ylabel 'u (m/s^2)';
set(f, 'Position', [2200 800 560 180]);
axis([0 8 -5 2.5]);
legend({'OPCBF', 'ECBF', 'MPC'},'Location','SouthEast','FontSize',12);

%% Constraint Plots
f = figure(4); clf;
plot(data1.t(1:end-1), data1.h, 'b', 'LineWidth', 1.5); hold on;
plot(data2.t(1:end-1), data2.h, 'g', 'LineWidth', 1.5);
plot(data3.t(1:end-1), data3.h, 'r', 'LineWidth', 1.5);
xlabel 'Time (s)';
ylabel 'h (meters)';
set(f, 'Position', [2200 500 560 180]);
axis([2.4 4.4 -20 0]);
legend({'OPCBF', 'ECBF', 'MPC'},'Location','SouthEast','FontSize',12);

f = figure(5); clf;
plot(data4.t(1:end-1), data4.h, 'b', 'LineWidth', 1.5); hold on;
plot(data5.t(1:end-1), data5.h, 'g', 'LineWidth', 1.5);
plot(data6.t(1:end-1), data6.h, 'r', 'LineWidth', 1.5);
xlabel 'Time (s)';
ylabel 'h (meters)';
set(f, 'Position', [2200 200 560 180]);
axis([2.4 4.4 -20 0]);
legend({'OPCBF', 'ECBF', 'MPC'},'Location','SouthEast','FontSize',12);

%% CBF Plots
f = figure(6); clf;
plot(data1.t(1:end-1), data1.H, 'b', 'LineWidth', 1.5); hold on;
plot(data2.t(1:end-1), data2.H, 'g', 'LineWidth', 1.5);
plot(data3.t(1:end-1), data3.h, 'r', 'LineWidth', 1.5);
xlabel 'Time (s)';
ylabel 'H';
set(f, 'Position', [2800 500 560 180]);
axis([2.4 4.4 -12 0]);
legend({'OPCBF', 'ECBF', 'MPC'},'Location','SouthEast','FontSize',12);

f = figure(7); clf;
plot(data4.t(1:end-1), data4.H, 'b', 'LineWidth', 1.5); hold on;
plot(data5.t(1:end-1), data5.H, 'g', 'LineWidth', 1.5);
plot(data6.t(1:end-1), data6.h, 'r', 'LineWidth', 1.5);
xlabel 'Time (s)';
ylabel 'H';
set(f, 'Position', [2800 200 560 180]);
axis([2.4 4.4 -12 0]);
legend({'OPCBF', 'ECBF', 'MPC'},'Location','SouthEast','FontSize',12);

%% Videos
global sim_case
sim_case = 1;

record = 1;
if record
    vidObj = VideoWriter('movie.avi');
    vidObj.FrameRate = 20; % double real speed
    open(vidObj);
end

f = figure(1); clf;
l1 = lane1(data2.x(1,1));
l2 = lane2(data2.x(3,1));
p1 = plot(l1(1), l1(2), 'bo', 'MarkerFaceColor', 'b'); hold on;
p2 = plot(l2(1), l2(2), 'go', 'MarkerFaceColor', 'g');
plot([-50, -3], [3, 3], 'k')
plot([-50, -3], [-3, -3], 'k')
plot([50, 3], [3, 3], 'k')
plot([50, 3], [-3, -3], 'k')
plot([-3, -3], [-50, -3], 'k')
plot([3, 3], [-50, -3], 'k')
plot([-3, -3], [50, 3], 'k')
plot([3, 3], [50, 3], 'k')
axis equal; axis([-50 50 -50 50]);
xlabel 'x (meters)';
ylabel 'y (meters)';
title 'Straight Intersection - Existing Approach (ECBF)';
set(f, 'Position', [1800 200 520 420]);
set(gca, 'XTick', -50:10:50);
if record
    for i=1:6
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end
N = length(data2.x);
for i=1:10:N
    l1 = lane1(data2.x(1,i));
    l2 = lane2(data2.x(3,i));
    set(p1, 'XData', l1(1), 'YData', l1(2));
    set(p2, 'XData', l2(1), 'YData', l2(2));
    if record
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end

l1 = lane1(data1.x(1,1));
l2 = lane2(data1.x(3,1));
set(p1, 'XData', l1(1), 'YData', l1(2));
set(p2, 'XData', l2(1), 'YData', l2(2));
title 'Straight Intersection - New Approach (OPCBF)';
if record
    for i=1:6
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end
N = length(data1.x);
for i=1:10:N
    l1 = lane1(data1.x(1,i));
    l2 = lane2(data1.x(3,i));
    set(p1, 'XData', l1(1), 'YData', l1(2));
    set(p2, 'XData', l2(1), 'YData', l2(2));
    if record
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end

l1 = lane1(data3.x(1,1));
l2 = lane2(data3.x(3,1));
set(p1, 'XData', l1(1), 'YData', l1(2));
set(p2, 'XData', l2(1), 'YData', l2(2));
title 'Straight Intersection - Comparison to MPC';
if record
    for i=1:6
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end
N = length(data3.x);
for i=1:N
    l1 = lane1(data3.x(1,i));
    l2 = lane2(data3.x(3,i));
    set(p1, 'XData', l1(1), 'YData', l1(2));
    set(p2, 'XData', l2(1), 'YData', l2(2));
    if record
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end

sim_case = 2;
l1 = lane1(data5.x(1,1));
l2 = lane2(data5.x(3,1));
set(p1, 'XData', l1(1), 'YData', l1(2));
set(p2, 'XData', l2(1), 'YData', l2(2));
title 'Left Turn - Existing Approach (ECBF)';
if record
    for i=1:6
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end
N = length(data5.x);
for i=1:10:N
    l1 = lane1(data5.x(1,i));
    l2 = lane2(data5.x(3,i));
    set(p1, 'XData', l1(1), 'YData', l1(2));
    set(p2, 'XData', l2(1), 'YData', l2(2));
    if record
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end

l1 = lane1(data4.x(1,1));
l2 = lane2(data4.x(3,1));
set(p1, 'XData', l1(1), 'YData', l1(2));
set(p2, 'XData', l2(1), 'YData', l2(2));
title 'Left Turn - New Approach (OPCBF)';
if record
    for i=1:6
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end
N = length(data4.x);
for i=1:10:N
    l1 = lane1(data4.x(1,i));
    l2 = lane2(data4.x(3,i));
    set(p1, 'XData', l1(1), 'YData', l1(2));
    set(p2, 'XData', l2(1), 'YData', l2(2));
    if record
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end

l1 = lane1(data6.x(1,1));
l2 = lane2(data6.x(3,1));
set(p1, 'XData', l1(1), 'YData', l1(2));
set(p2, 'XData', l2(1), 'YData', l2(2));
title 'Left Turn - Comparison to MPC';
if record
    for i=1:6
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end
N = length(data6.x);
for i=1:N
    l1 = lane1(data6.x(1,i));
    l2 = lane2(data6.x(3,i));
    set(p1, 'XData', l1(1), 'YData', l1(2));
    set(p2, 'XData', l2(1), 'YData', l2(2));
    if record
        frame = getframe(f);
        writeVideo(vidObj,frame);
    end
end

if record
    close(vidObj);
end