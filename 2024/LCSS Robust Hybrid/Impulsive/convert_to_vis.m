function [x_hcw, rscale, vscale, v_vec] = convert_to_vis(t, x_eci)
% Converts a 4xN vector of inertial frame coordinates to a fake-hcw frame. The obstacle positions
% are *nearly* constant in the hcw frame, so this is helpful for visualization of trajectories.
mu = 398600e9;
[c, ~, oe] = get_center(0);
n = sqrt(mu/oe.a^3);
rscale = 1./vecnorm(x_eci(1:2,:));
vscale = 1./sqrt(mu./rscale);
v_vec = c(3:4)/norm(c(3:4));
x_hcw = zeros(size(x_eci));
for i=1:length(t)
    [xc, theta] = get_center(t(i));
    R = [cos(theta), sin(theta); -sin(theta), cos(theta)];
    x_hcw(1:2,i) = R*(x_eci(1:2,i) - xc(1:2));
%     x_hcw(1:2,i) = x_hcw(1:2,i)*rscale(i); % standardized position reference
    
    if size(x_eci, 1) == 4
        x_hcw(3:4,i) = R*(x_eci(3:4,i) - xc(3:4)) - [eye(2), zeros(2,1)]*cross([0;0;n], [x_hcw(1:2,i); 0]);
%         x_hcw(3:4,i) = x_hcw(3:4,i)*vscale(i); % standardized velocity reference
    end
end
end