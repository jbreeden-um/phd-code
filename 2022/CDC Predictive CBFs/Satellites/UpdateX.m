function xr = UpdateX(t1, t2, x1, u)
global time_step
[~, xr] = ode45(@(t,x) two_body(x,u), t1:time_step:t2, x1, odeset('RelTol', 1e-8, 'AbsTol', 1e-8));
end