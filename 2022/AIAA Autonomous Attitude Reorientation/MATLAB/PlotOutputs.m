function stop = PlotOutputs(t, x, extra)
global s_target ctheta
persistent PositionPlotHandle1 CurrentPositionHandle1 constants ...
    PositionPlotHandle2 CurrentPositionHandle2 ExclusionPlotHandle

if t(1)==0
    constants = load('parameters.mat');
end

if isequal(extra, 'init')
    figure(1); clf;
    ExclusionPlotHandle = PlotCone(get_s(t(1)), acos(ctheta));
    hold on;
    [s1, s2, s3] = sphere(20);
    surf(s1, s2, s3, 'FaceAlpha', 0);
    
    p1 = QtoR(x(1,1:4))*constants.p1;
    p2 = QtoR(x(1,1:4))*constants.p2;
    PositionPlotHandle1 = plot3(p1(1), p1(2), p1(3), 'b', 'LineWidth', 3);
    CurrentPositionHandle1 = plot3(p1(1), p1(2), p1(3), 'bo', 'MarkerFaceColor', [0;0;1]);
    PositionPlotHandle2 = plot3(p2(1), p2(2), p2(3), 'c', 'LineWidth', 3);
    CurrentPositionHandle2 = plot3(p2(1), p2(2), p2(3), 'co', 'MarkerFaceColor', [0;1;1]);
    
    TargetPositionHandle = plot3(s_target(1), s_target(2), s_target(3), 'go', 'MarkerFaceColor', [0;1;0]);
    xlabel 'x'; ylabel 'y'; zlabel 'z'; hold on;
    axis equal;
elseif isequal(extra, 'done')
    % Figure out if you need any cleanup here
else
    PlotCone(get_s(t(end)), acos(ctheta), ExclusionPlotHandle);
    p1 = RotQ(constants.p1, x(:,1:4)');
    p2 = RotQ(constants.p2, x(:,1:4)');
    set(PositionPlotHandle1, 'XData', p1(1,:), 'YData', p1(2,:), 'ZData', p1(3,:));
    set(CurrentPositionHandle1, 'XData', p1(1,end), 'YData', p1(2,end), 'ZData', p1(3,end));
    set(PositionPlotHandle2, 'XData', p2(1,:), 'YData', p2(2,:), 'ZData', p2(3,:));
    set(CurrentPositionHandle2, 'XData', p2(1,end), 'YData', p2(2,end), 'ZData', p2(3,end));
end

stop = 0;
end