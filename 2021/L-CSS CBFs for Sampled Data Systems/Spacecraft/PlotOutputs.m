function stop = PlotOutputs(t, x, extra)
global outdata s_target separation two_obstacles constants
persistent PositionPlotHandle CurrentPositionHandle TargetPositionHandle

if isequal(extra, 'init')
    figure(1); clf;
    PlotCone(constants.s, constants.theta);
    hold on;
    [s1, s2, s3] = sphere(20);
    surf(s1, s2, s3, 'FaceAlpha', 0);
    
    PositionPlotHandle = plot3(x(:,1), x(:,2), x(:,3), 'b', 'LineWidth', 3);
    CurrentPositionHandle = plot3(x(end,1), x(end,2), x(end,3), 'bo', 'MarkerFaceColor', [0;0;1]);
    TargetPositionHandle = plot3(s_target(1), s_target(2), s_target(3), 'go', 'MarkerFaceColor', [0;1;0]);
    xlabel 'x'; ylabel 'y'; zlabel 'z'; hold on;
    axis equal;
elseif isequal(extra, 'done')
    % Figure out if you need any cleanup here
else
    set(PositionPlotHandle, 'XData', x(:,1), 'YData', x(:,2), 'ZData', x(:,3));
    set(CurrentPositionHandle, 'XData', x(end,1), 'YData', x(end,2), 'ZData', x(end,3));
end

stop = 0;
end