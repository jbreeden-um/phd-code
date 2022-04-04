function xdot = two_body(x,u)
r = x(1:3);
v = x(4:6);
mu = 398600;
xdot = [v; -mu*r/norm(r)^3 + u];
end