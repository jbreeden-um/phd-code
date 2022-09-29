%% Initial Problem
rb = [cosd(45); sind(45); 0];
r1 = [1; 0; 0];
r2 = [0; 1; 0];
theta = 50;
M1 = [r1*rb' + rb*r1' - (rb'*r1 + cosd(theta))*eye(3), cross(rb, r1); cross(rb, r1)', rb'*r1 - cosd(theta)];
M2 = [r2*rb' + rb*r2' - (rb'*r2 + cosd(theta))*eye(3), cross(rb, r2); cross(rb, r2)', rb'*r2 - cosd(theta)];

figure(1); clf;
[s1, s2, s3] = sphere(20);
surf(s1, s2, s3, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'EdgeAlpha', 0.5); hold on; axis equal;
[c1, c2, c3] = cylinder(sind(theta), 40);
c1(1,:) = 0;
c2(1,:) = 0;
c3(2,:) = c3(2,:)*cosd(theta);
c_all = [c1(:)'; c2(:)'; c3(:)'];
c_new = RotQ(c_all, [0; sin(pi/4); 0; cos(pi/4)]); % rotate the +z axis to the +x axis
d1 = reshape(c_new(1,:), size(c1));
d2 = reshape(c_new(2,:), size(c2));
d3 = reshape(c_new(3,:), size(c3));
surf(d1, d2, d3, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'LineWidth', 2, 'EdgeColor', 'r', 'EdgeAlpha', 0);

c_new = RotQ(c_all, [-sin(pi/4); 0; 0; cos(pi/4)]); % rotate the +z axis to the +y axis
d1 = reshape(c_new(1,:), size(c1));
d2 = reshape(c_new(2,:), size(c2));
d3 = reshape(c_new(3,:), size(c3));
surf(d1, d2, d3, 'FaceColor', [1; 0.5; 0], 'FaceAlpha', 0.5, 'LineWidth', 2, 'EdgeColor', 'r', 'EdgeAlpha', 0);

xlabel x; ylabel y;
view(132, 16);

%% Expert Solution
q0 = [0.18; -0.18; 0; 0.96];
q_crit = fsolve(@(q) [q'*M1*q; q'*M2*q; norm(q)-1], q0);

vc = RotQ(rb, q_crit);
% plot3([0 vc(1)], [0 vc(2)], [0 vc(3)], 'm', 'LineWidth', 2);
theta3 = round(asind(vc(3)))+1;
r3 = [vc(1:2); 0]; r3 = r3/norm(r3);
M3 = [r3*rb' + rb*r3' - (rb'*r3 + cosd(theta3))*eye(3), cross(rb, r3); cross(rb, r3)', rb'*r3 - cosd(theta3)];

if 0
[c1, c2, c3] = cylinder(sind(theta3), 40);
c1(1,:) = 0;
c2(1,:) = 0;
c3(2,:) = c3(2,:)*cosd(theta3);
c_all = [c1(:)'; c2(:)'; c3(:)'];
c_new = RotQ(c_all, [-sin(pi/4)/sqrt(2); sin(pi/4)/sqrt(2); 0; cos(pi/4)]); % rotate the +z axis to the +x axis
d1 = reshape(c_new(1,:), size(c1));
d2 = reshape(c_new(2,:), size(c2));
d3 = reshape(c_new(3,:), size(c3));
surf(d1, d2, d3, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'LineWidth', 2, 'EdgeColor', 'r', 'EdgeAlpha', 0);
end

%% Numerical Solution
% vset = [0.636100 0.618349 -0.461543 0.034907;
%         0.636100 0.618349 0.461543 0.034907]; % Case 1
% vset = [0.707815  0.704394  -0.053169 0.436332;
%         0.711703  0.700292  0.055401 0.436332]; % Case 2
% vset = [0.706744  0.707451  0.005064 0.471239]; % Case 3
vset = load('data/vectors.dat');

for i=1:size(vset, 1)
    v = vset(i,1:3)';
    theta3 = vset(i,4)*180/pi;
    M3 = [v*rb' + rb*v' - (rb'*v + cosd(theta3))*eye(3), cross(rb, v); cross(rb, v)', rb'*v - cosd(theta3)];
    [c1, c2, c3] = cylinder(sind(theta3), 40);
    c1(1,:) = 0;
    c2(1,:) = 0;
    c3(2,:) = c3(2,:)*cosd(theta3);
    c_all = [c1(:)'; c2(:)'; c3(:)'];
    a = cross([0;0;1],v);
    phi = acos(v(3)/norm(v));
    c_new = RotQ(c_all, [a/norm(a)*sin(phi/2); cos(phi/2)]);
    d1 = reshape(c_new(1,:), size(c1));
    d2 = reshape(c_new(2,:), size(c2));
    d3 = reshape(c_new(3,:), size(c3));
    surf(d1, d2, d3, 'FaceColor', [.5; 0.25; 0], 'FaceAlpha', 0.5, 'LineWidth', 2, 'EdgeColor', 'r', 'EdgeAlpha', 0);
end

colors = 'bg';
for i=1:size(vset,1)
    cluster = load(['data/cluster' num2str(i+1) '.dat']);
    vcluster = zeros(3,length(cluster));
    for j=1:length(cluster)
        vcluster(:,j) = RotQ(rb, cluster(j,1:4)');
    end
    plot3(vcluster(1,:), vcluster(2,:), vcluster(3,:), '.', 'Color', colors(i));
end