function xdot = two_body(t,x,u)
if nargin==2, u=[0;0;0]; end
mu = 398600e9;
xdot = [x(4:6); -mu*x(1:3)/norm(x(1:3))^3+u(1:3)];
end