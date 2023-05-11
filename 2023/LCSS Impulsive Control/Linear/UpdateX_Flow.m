function [t,x] = UpdateX_Flow(t0, x0, T, N_to_return)
global A_sys
A = A_sys;
t_span = linspace(t0, t0+T, N_to_return);
[t_sim,x_sim] = ode45(@(t,x) A*x, t_span, x0);
t = t_sim';
x = x_sim';
end