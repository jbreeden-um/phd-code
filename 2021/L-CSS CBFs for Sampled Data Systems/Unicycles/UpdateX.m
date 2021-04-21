function x1 = UpdateX(t, x0, u, dt)
% Implements the dynamics in a ZOH fashion

x1 = zeros(3,1);
if abs(u(2)) > 1e-10
    x1(1) = x0(1) + u(1)/u(2)*(sin(x0(3) + u(2)*dt) - sin(x0(3)));
    x1(2) = x0(2) + u(1)/u(2)*(-cos(x0(3) + u(2)*dt) + cos(x0(3)));
else
    x1(1) = x0(1) + u(1)*cos(x0(3))*dt;
    x1(2) = x0(2) + u(1)*sin(x0(3))*dt;
end

x1(3) = wrap_angle(x0(3) + u(2)*dt);
end