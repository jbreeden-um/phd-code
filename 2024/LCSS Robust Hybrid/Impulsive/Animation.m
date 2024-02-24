% Animation that you can run right after a simulation to visualize results

figure(11); clf;
theta = linspace(0, 2*pi, 100);
cx = cos(theta);
cy = sin(theta);
rho = 1;
xc = get_center(0);
obs = obstacle_location(1,0); o1 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r'); hold on;
obs = obstacle_location(2,0); o2 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(3,0); o3 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
obs = obstacle_location(4,0); o4 = fill(obs(1)/1e3+cx*rho, obs(2)/1e3+cy*rho, 'r');
dock_width = 100;
o5 = fill(xc(1)/1e3 + dock_width*[-1, 1, 1, -1, -1], xc(2)/1e3 + dock_width*[0, 0, 1, 1, 0], 'r');
o_array = [o1; o2; o3; o4];
xlabel 'x_1'; ylabel 'x_2';
axis([xc(1)/1e3-10, xc(1)/1e3+10, xc(2)/1e3-10, xc(2)/1e3+10]);
axis equal;
px = plot(x(1,1), x(2,1), 'ko', 'MarkerFaceColor', 'b');
pc = plot(xc(1), xc(2), 'k.', 'MarkerSize', 6);

for i=1:length(t)
    xi = x(:,i);
    for j=1:4
        obs = obstacle_location(j,t(i));
        set(o_array(j), 'Vertices', [obs(1)/1e3+cx(:)*rho, obs(2)/1e3+cy(:)*rho]);
    end
    [xc, theta] = get_center(t(i));
    rectangle = [-1, 1, 1, -1, -1; 0, 0, 1, 1, 0];
    rectangle = [cos(theta), -sin(theta); sin(theta), cos(theta)]*rectangle;
    set(o5, 'Vertices', [xc(1), xc(2)]/1e3 + dock_width*rectangle');
    set(px, 'XData', xi(1)/1e3, 'YData', xi(2)/1e3);
    set(pc, 'XData', xc(1)/1e3, 'YData', xc(2)/1e3);
    axis([xc(1)/1e3-20, xc(1)/1e3+20, xc(2)/1e3-20, xc(2)/1e3+20]);
    drawnow;
%     pause(0.01);
end