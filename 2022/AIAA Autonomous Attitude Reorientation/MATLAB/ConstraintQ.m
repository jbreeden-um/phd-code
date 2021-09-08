function [A, b] = ConstraintQ(t, x, constants, p, index)
% Returns the A and b constraint matrices for ensuring avoidance is achieved.
% Results should be fed into QP as A*u \leq b
global outdata

outdata.h(index,1) = h(t,x,p);
outdata.hdot(index,1) = hdot(t,x,p);

H = outdata.h(index,1) + absSq(outdata.hdot(index,1))/(2*constants.mu);
if H > 0
    disp 'Excessive H'
end
if outdata.h > 0
    disp 'Excessive h';
end
outdata.H(index,1) = H;

dt = constants.dt;
M2plus = constants.M2plus;
M3plus = constants.M3plus;
h_future = @(phi_val) h(t, x, p) + hdot(t, x, p)*dt + 1/2*phi_val*dt^2 + 1/2*M2plus*dt^2 + 1/6*M3plus*dt^3;
phi_req1 = fzero(@(y) h_future(y) + constants.delta_2, 0);
b1 = phi_req1 - phi(t,x,[0;0;0;0],p);

hdot_future = @(phi_val) hdot(t, x, p) + phi_val*dt + M2plus*dt + 1/2*M3plus*dt^2;
H_future = @(phi_val) h_future(phi_val) + 1/(2*constants.mu)*absSq(hdot_future(phi_val));
phi_req2 = fzero(@(y) H_future(y) + constants.Delta_2, 0);
b2 = phi_req2 - phi(t,x,[0;0;0;0],p);

b = min(b1, b2);
A = [0,0,0,0];
for i=1:4
    dist = 1; % arbitrary since it is linear
    e = [0;0;0;0];
    f = e; f(i) = dist;
    A(i) = (phi(t,x,f,p) - phi(t,x,e,p))/dist;
end
% Using a numerical gradient for phi because it is linear rather than writing another 
% function to compute it.
end