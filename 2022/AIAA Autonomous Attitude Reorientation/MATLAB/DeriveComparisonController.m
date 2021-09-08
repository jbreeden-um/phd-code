% This file derives the formulas for the gradients and hessians required to implement the
% control law being used for comparison.

syms q0 q1 q2 q3 real
syms qr0 qr1 qr2 qr3 real
q = [q0; q1; q2; q3];
qr = [qr0; qr1; qr2; qr3];
term1 = QxQ(QConj(qr), q) - [1;0;0;0];
term1 = term1.'*term1;

% grad = jacobian(term1, q)

syms k
syms x11 x12 x13 real
syms y11 y12 y13 y21 y22 y23 real
syms ctheta1 ctheta2
x1 = [x11; x12; x13];
y1 = [y11; y12; y13];
y2 = [y21; y22; y23];
M1 = [dot(x1, y1) - ctheta1, cross(x1, y1).';
    cross(x1, y1), x1*y1.' + y1*x1.' - (dot(x1,y1)+ctheta1)*eye(3)];
M2 = [dot(x1, y2) - ctheta2, cross(x1, y2)';
    cross(x1, y2), x1*y2.' + y2*x1.' - (dot(x1,y2)+ctheta2)*eye(3)];

Q = -eye(4); Q(1,1) = 1;
M1 = Q*M1*Q; M2 = Q*M2*Q;

V = term1*(-k*log(-q.'*M1*q/2) - k*log(-q.'*M2*q/2));

grad = jacobian(V, q).'
hess = hessian(V, q)