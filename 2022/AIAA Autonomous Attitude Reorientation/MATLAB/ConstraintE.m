function [H, A, b] = ConstraintE(t, x, constants, index)
% Returns the H, A, b constraint matrices for ensuring energy constraint is achieved.
% Results should be fed into QCQP as u'*H*u + A*u + b \leq 0
global outdata

h_omega = @(w) w'*constants.Je*w - constants.E_max;
w = x(5:7);
outdata.h(index,1) = h_omega(w);

if outdata.h(index,1) > 0
    disp 'Excessive kinetic energy';
end

% Quadratic Case
Z = constants.Z;
dt = constants.dt;
M2omega = constants.M2omega;
M1omega = constants.M1omega;
H = 1/2*2*Z(1:3,4:7)'*constants.Je*Z(1:3,4:7)*dt^2;
A = 2*w'*constants.Je*Z(1:3,4:7)*dt;
b = -h_omega(w) - M1omega*dt - 1/2*M2omega*dt^2;

% Linear Case
% M2omega_lin = constants.M2omega_lin;
% H = zeros(4);
% b = -h_omega(w) - M1omega*dt - 1/2*M2omega_lin*dt^2;
end