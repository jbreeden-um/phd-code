function x_hcw = convert_to_hcw(t, x_eci, frame_only)
if nargin==2
    frame_only = 0;
end
% Converts a 6xN vector of inertial frame coordinates to hcw frame.
mu = 398600e9;
[~, oe] = get_center(0);
n = sqrt(mu/oe.a^3);
x_hcw = zeros(size(x_eci));
R1 = [cos(oe.Omega), -sin(oe.Omega), 0; sin(oe.Omega), cos(oe.Omega), 0; 0, 0, 1]...
     *[1, 0, 0; 0, cos(oe.i), -sin(oe.i); 0, sin(oe.i), cos(oe.i)];
for i=1:length(t)
    [xc, oe] = get_center(t(i));
    theta = oe.omega + oe.nu;
    RP2I = R1*[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    RI2P = RP2I';
    if frame_only
        x_hcw(1:3,i) = RI2P*x_eci(1:3,i);
    else
        x_hcw(1:3,i) = RI2P*(x_eci(1:3,i) - xc(1:3));
        x_hcw(4:6,i) = RI2P*(x_eci(4:6,i) - xc(4:6)) - cross([0;0;n], x_hcw(1:3,i));
    end
end
end