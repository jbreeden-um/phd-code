function x_hcw = convert_to_hcw(t, x_eci)
% Converts a 4xN vector of inertial frame coordinates to hcw frame. The obstacle positions
% are *nearly* constant in the hcw frame, so this is helpful for visualization of
% trajectories.
mu = 398600e9;
Re = 6378e3;
alt = 600e3;
n = sqrt(mu/(Re+alt)^3);
x_hcw = zeros(size(x_eci));
for i=1:length(t)
    [xc, theta] = get_center(t(i));
    R = [cos(theta), sin(theta); -sin(theta), cos(theta)];
    x_hcw(1:2,i) = R*(x_eci(1:2,i) - xc(1:2));
    x_hcw(3:4,i) = R*(x_eci(3:4,i) - xc(3:4)) - [eye(2), zeros(2,1)]*cross([0;0;n], [x_hcw(1:2,i); 0]);
end
end