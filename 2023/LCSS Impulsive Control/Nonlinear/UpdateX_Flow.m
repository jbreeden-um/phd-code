function [t,x] = UpdateX_Flow(t0, x0, T, N_to_return)
t_span = linspace(t0, t0+T, N_to_return);
[t_sim,x_sim] = ode45(@two_body, t_span, x0);
t = t_sim';
x = x_sim';
end