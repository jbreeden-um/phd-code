% This file simulates the obstacle locations over time in an inertial frame, so that we
% can see what their orbits correspond to

mu = 398600e9;
rho = 1e3;
theta = linspace(0, 2*pi, 100);
cx = rho*cos(theta);
cy = rho*sin(theta);

[x0, ~, oe] = get_center(0);

figure(1); clf;
p_center = plot(x0(1), x0(2), 'ko', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
nu = linspace(0, 2*pi, 10000);
ox = zeros(2, length(nu));
for j=1:length(nu)
    [or, ov] = KeplerToCartesian(oe.a, oe.e, oe.i, oe.Omega, oe.omega, nu(j), mu);
    ox(:,j) = or(1:2);
end
plot(ox(1,:), ox(2,:), 'b', 'LineWidth', 2);

num_obs = 6;
for i=1:num_obs
    [x0, oe] = obstacle_location(i,0);
    nu = linspace(0, 2*pi, 10000);
    ox = zeros(2, length(nu));
    for j=1:length(nu)
        [or, ov] = KeplerToCartesian(oe.a, oe.e, oe.i, oe.Omega, oe.omega, nu(j), mu);
        ox(:,j) = or(1:2);
    end
    plot(ox(1,:), ox(2,:), 'k');
    p_obs{i} = fill(x0(1) + cx, x0(2) + cy, 'r');
end
axis equal;

period = 2*pi*sqrt(oe.a^3/mu);
tspan = 0:10:ceil(period);
% tspan = 0:100:00;
for t=tspan
    xc = get_center(t);
    bds = [xc(1), xc(1), xc(2), xc(2)];
    set(p_center, 'XData', xc(1), 'YData', xc(2));
    for i=1:num_obs
        x = obstacle_location(i,t);
        set(p_obs{i}, 'Vertices', [cx(:), cy(:)] + [x(1) x(2)])
        
        bds(1) = min(bds(1), x(1));
        bds(2) = max(bds(2), x(1));
        bds(3) = min(bds(3), x(2));
        bds(4) = max(bds(4), x(2));
    end
    bds = bds + [-2e3, 2e3, -2e3, 2e3];
    axis(bds);
    
    drawnow;
    
    if max(abs(bds - [xc(1), xc(1), xc(2), xc(2)])) > 100e3
        something_is_wrong = 1;
    end
end