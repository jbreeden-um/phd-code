function x1 = UpdateX(t, x0, u, dt)
% Implements the dynamics in a ZOH fashion

x1 = zeros(6,1);
x1(1:3) = expm(skew(x0(4:6)*dt + 1/2*u*dt^2))*x0(1:3);
x1(4:6) = x0(4:6) + u*dt;
end