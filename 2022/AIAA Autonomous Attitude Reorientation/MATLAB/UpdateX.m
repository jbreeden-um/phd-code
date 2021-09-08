function x1 = UpdateX(t, x0, u, dt)
% Implements the dynamics in a ZOH fashion
persistent constants v0;
if t==0
    constants = load('parameters.mat');
    rngstate = load('rngstate.mat');
    rng(rngstate.rngstate);
    v0 = rand(3,1)*constants.torque;
    if norm(v0) > constants.torque, v0 = v0/norm(v0)*constants.torque; end
end

% q = x(1:4);
% w = x(5:7);
% W = x(8:11);
Z = constants.Z;
Jtot = constants.Jtot;
Jw = constants.Jw;
Aw = constants.wheel_axis;
Q = @(w) [0, -w(1), -w(2), -w(3);
    w(1), 0, w(3), -w(2);
    w(2), -w(3), 0, w(1);
    w(3), w(2), -w(1), 0];
v1 = rand(3,1)*constants.torque; if norm(v1) > constants.torque, v1 = v1/norm(v1)*constants.torque; end
v = @(tau) v0 + (v1-v0)*(tau - t); % disturbance function
    % Disturbance function is regulated to be continuous.
f = @(t, x) [1/2*Q(x(5:7))*x(1:4);
    Z*[v(t) - cross(x(5:7), Jtot*x(5:7) + Aw*Jw*x(8:11)); u]];
[~,y] = ode45(f, [t t+dt], x0);
x1 = y(end,:)';

v0 = v1;
end