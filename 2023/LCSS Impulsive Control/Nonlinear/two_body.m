function xdot = two_body(t,x)
mu = 398600e9;
xdot = [x(3:4); -mu*x(1:2)/norm(x(1:2))^3];
end