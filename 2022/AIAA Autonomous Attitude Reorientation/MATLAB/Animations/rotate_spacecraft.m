function spacecraft = rotate_spacecraft(R, spacecraft)
mesh = spacecraft.meshes;
surf = spacecraft.surfaces;
line = spacecraft.lines;
for i=1:length(mesh)
    set(mesh(i), 'Vertices', mesh(i).Vertices*R');
end
for i=1:length(surf)
    % Only works if x data is a matrix rather than a vector
    sx = size(surf(i).XData);
    x0 = reshape(surf(i).XData, [prod(sx), 1]);
    y0 = reshape(surf(i).YData, [prod(sx), 1]);
    z0 = reshape(surf(i).ZData, [prod(sx), 1]);
    w = [x0, y0, z0]*R';
    x1 = reshape(w(:,1), sx);
    x2 = reshape(w(:,2), sx);
    x3 = reshape(w(:,3), sx);
    set(surf(i), 'XData', x1);
    set(surf(i), 'YData', x2);
    set(surf(i), 'ZData', x3);
end
for i=1:length(line)
    x0 = line(i).XData(:);
    y0 = line(i).YData(:);
    z0 = line(i).ZData(:);
    w = [x0, y0, z0]*R';
    set(line(i), 'XData', w(:,1));
    set(line(i), 'YData', w(:,2));
    set(line(i), 'ZData', w(:,3));
end

pos = R*spacecraft.fixed_position;
delta = pos - spacecraft.fixed_position;
spacecraft.fixed_position = pos;
line = spacecraft.fixed_lines;
for i=1:length(line)
    set(line(i), 'XData', line(i).XData + delta(1));
    set(line(i), 'YData', line(i).YData + delta(2));
    set(line(i), 'ZData', line(i).ZData + delta(3));
end
surf = spacecraft.fixed_surfaces;
for i=1:length(surf)
    set(surf(i), 'XData', surf(i).XData + delta(1));
    set(surf(i), 'YData', surf(i).YData + delta(2));
    set(surf(i), 'ZData', surf(i).ZData + delta(3));
end

spacecraft.R = spacecraft.R*R;
end