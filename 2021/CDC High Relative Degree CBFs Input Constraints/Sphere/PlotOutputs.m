function stop = PlotOutputs(t, x, extra)
global outdata plates_to_consider ground_vec
persistent PositionPlotHandle GroundPlotHandle ClosestPointPlotHandle tStart tEnd Eros CurrentPositionHandle

if isequal(extra, 'init')
    tStart = t(1);
    tEnd = t(2);
    
    figure(1); clf;
    r = 10;
    [s1, s2, s3] = sphere(20);
    surf(s1*r, s2*r, s3*r, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeAlpha', 0.7);
    hold on;
    ClosestPointPlotHandle = plot3(nan, nan, nan, 'ro', 'MarkerFaceColor', [1;0;0]);
    PositionPlotHandle = plot3(x(:,1), x(:,2), x(:,3), 'b');
    CurrentPositionHandle = plot3(x(:,1), x(:,2), x(:,3), 'bo', 'MarkerFaceColor', [0;0;1]);
    GroundPlotHandle = plot3(nan, nan, nan, 'go', 'MarkerFaceColor', [0;1;0]);
    xlabel 'x'; ylabel 'y'; zlabel 'z'; hold on;
    axis equal;
elseif isequal(extra, 'done')
    % Figure out if you need any cleanup here
else
    set(PositionPlotHandle, 'XData', x(:,1), 'YData', x(:,2), 'ZData', x(:,3));
    set(CurrentPositionHandle, 'XData', x(end,1), 'YData', x(end,2), 'ZData', x(end,3));
    set(ClosestPointPlotHandle, 'XData', outdata.rc(1), 'YData', outdata.rc(2), 'ZData', outdata.rc(3));
    set(GroundPlotHandle, 'XData', ground_vec(1), 'YData', ground_vec(2), 'ZData', ground_vec(3));
end


stop = 0;
end