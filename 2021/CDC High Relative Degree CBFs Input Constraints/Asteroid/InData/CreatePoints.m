load Eros_Shape.mat

%% Examine Spline
np = length(plates);
nv = length(vertices);
distances = zeros(np, 3);

for i=1:np
    row = plates(i,:);
    v1 = vertices(row(1)+1,:);
    v2 = vertices(row(2)+1,:);
    v3 = vertices(row(3)+1,:);
    distances(i,1) = norm(v2-v1);
    distances(i,2) = norm(v3-v2);
    distances(i,3) = norm(v1-v3);
end

% Minimum distance between vertices in meters
min_dist = min(distances(:))*1e3

%% Make New Model
points = zeros(np, 3);
normals = zeros(np, 3);

for i=1:np
    row = plates(i,:);
    v1 = vertices(row(1)+1,:);
    v2 = vertices(row(2)+1,:);
    v3 = vertices(row(3)+1,:);
    points(i, :) = (v1 + v2 + v3)/3;
    normals(i, :) = cross(v3-v1, v2-v1);
    normals(i, :) = normals(i,:)/norm(normals(i,:))*sign(dot(points(i,:), normals(i,:)));
end

figure(1); clf;
trisurf(plates+1, vertices(:,1), vertices(:,2), vertices(:,3), ...
    'FaceColor', [.7; .7; .7], 'FaceAlpha', 0.4, 'EdgeAlpha', 0.1);
hold on;
for i=1:np
    x0 = points(i,:);
    x1 = x0 - normals(i, :);
    plot3([x0(1), x1(1)], [x0(2), x1(2)], [x0(3), x1(3)], 'k');
end
axis equal;

save('Eros_Points', 'points', 'normals');