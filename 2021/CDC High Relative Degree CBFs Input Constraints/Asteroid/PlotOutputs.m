function stop = PlotOutputs(t, x, extra)
global outdata plates_to_consider ground_vec
persistent PositionPlotHandle GroundPlotHandle PlatesPlotHandle tStart tEnd Eros CurrentPositionHandle

if isequal(extra, 'init')
    tStart = t(1);
    tEnd = t(2);
    
    figure(1); clf;
    Eros = load('InData/Eros_Shape.mat');
    trisurf(Eros.plates+1, Eros.vertices(:,1), Eros.vertices(:,2), Eros.vertices(:,3), ...
        'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1, 'EdgeAlpha', 0.7);
    hold on;
    PlatesPlotHandle = trisurf([], Eros.vertices(:,1), Eros.vertices(:,2), Eros.vertices(:,3), ...
        'FaceColor', [1 0 0], 'FaceAlpha', 1, 'EdgeAlpha', 0);
    PositionPlotHandle = plot3(x(:,1), x(:,2), x(:,3), 'b', 'LineWidth', 2);
    CurrentPositionHandle = plot3(x(:,1), x(:,2), x(:,3), 'bo', 'MarkerFaceColor', [0;0;1]);
    GroundPlotHandle = plot3(nan, nan, nan, 'go', 'MarkerFaceColor', [0;1;0]);
    xlabel 'x'; ylabel 'y'; zlabel 'z'; hold on;
    axis equal;
elseif isequal(extra, 'done')
    % Figure out if you need any cleanup here
else
    set(PositionPlotHandle, 'XData', x(:,1), 'YData', x(:,2), 'ZData', x(:,3));
    set(CurrentPositionHandle, 'XData', x(end,1), 'YData', x(end,2), 'ZData', x(end,3));
    set(PlatesPlotHandle, 'Faces', Eros.plates(plates_to_consider,:)+1);
%     disp(sum(plates_to_consider))
    set(GroundPlotHandle, 'XData', ground_vec(1), 'YData', ground_vec(2), 'ZData', ground_vec(3));
end




stop = 0;
end