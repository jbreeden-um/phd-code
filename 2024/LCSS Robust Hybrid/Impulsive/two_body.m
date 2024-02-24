function xdot = two_body(t,x,u)
if nargin==2, u=[0;0]; end
mu = 398600e9;
xdot = [x(3:4); -mu*x(1:2)/norm(x(1:2))^3+u(1:2)];
end