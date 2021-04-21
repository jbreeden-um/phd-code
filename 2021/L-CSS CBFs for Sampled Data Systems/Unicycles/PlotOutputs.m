function stop = PlotOutputs(t, x, extra)
global outdata r_target constants separation two_obstacles
persistent PositionPlotHandle CurrentPositionHandle TargetPositionHandle ...
    agent_radius

if isequal(extra, 'init')
    figure(1); clf;
    agent_radius = 1;
    obstacle_radius = constants.rho - agent_radius;
    
    x = linspace(-obstacle_radius, obstacle_radius, 100);
    y = sqrt(obstacle_radius^2 - x.^2);
    fill([x, fliplr(x)], [y, -y], 'r', 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
    hold on;
    if two_obstacles
        fill([x, fliplr(x)]-separation/sqrt(2), [y, -y]+separation/sqrt(2), 'r', 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
    end
    
    PositionPlotHandle = plot(x(:,1), x(:,2), 'b');
    [rx, ry] = draw_unicycle(x(1,:), agent_radius);
    CurrentPositionHandle = plot(rx, ry, 'b');
    TargetPositionHandle = plot(r_target(1), r_target(2), 'go', 'MarkerFaceColor', [0;1;0]);
    xlabel 'x'; ylabel 'y'; zlabel 'z'; hold on;
    axis equal;
elseif isequal(extra, 'done')
    % Figure out if you need any cleanup here
else
    set(PositionPlotHandle, 'XData', x(:,1), 'YData', x(:,2));
    [rx, ry] = draw_unicycle(x(end,:), agent_radius);
    set(CurrentPositionHandle, 'XData', rx, 'YData', ry);
end

stop = 0;
end

function [rx, ry] = draw_unicycle(x, radius)
n = 20;
phi = x(3);
theta = linspace(0, 2*pi, n);
theta = phi + [theta, 5*pi/6, 0, 7*pi/6, 2*pi];

rx = radius*cos(theta); rx(end-2) = 0;
ry = radius*sin(theta); ry(end-2) = 0;

rx = x(1) + rx;
ry = x(2) + ry;
end