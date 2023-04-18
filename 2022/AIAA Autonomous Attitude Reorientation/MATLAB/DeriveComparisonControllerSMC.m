% This file derives the formulas for the gradients required to implement the
% control law being used for comparison in ComparisonControllerSMC

syms q0 q1 q2 q3 real
syms qr0 qr1 qr2 qr3 real
q = [q0; q1; q2; q3];
qr = [qr0; qr1; qr2; qr3];
term1 = (qr - q).'*(qr - q);

syms k alpha real
syms x11 x12 x13 real
syms y11 y12 y13 y21 y22 y23 real
syms ctheta1 ctheta2
x1 = [x11; x12; x13]; % inertial coordinates
y1 = [y11; y12; y13]; % body fixed coordinates
y2 = [y21; y22; y23];

% This paper suggests that this should be cross(y, x) instead, but that
% results in simulation error. It should be cross(x, y), as is done below
% and in the barrier functions paper.
M1 = [dot(x1, y1) - ctheta1, cross(x1, y1).';
    cross(x1, y1), x1*y1.' + y1*x1.' - (dot(x1,y1)+ctheta1)*eye(3)];
M2 = [dot(x1, y2) - ctheta2, cross(x1, y2)';
    cross(x1, y2), x1*y2.' + y2*x1.' - (dot(x1,y2)+ctheta2)*eye(3)];

Q = -eye(4); Q(1,1) = 1;
M1 = Q*M1*Q; M2 = Q*M2*Q;

V = term1*alpha*(1/(q.'*M1*q)^2 + 1/(q.'*M2*q)^2);

grad = jacobian(V, q).'