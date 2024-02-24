function [t,x] = UpdateX_Flow(t0, x0, T, N_to_return, use_disturbance)
if nargin==5 && use_disturbance
    t = linspace(t0, t0+T, N_to_return);
    x = my_ode4(t, x0);
else
    t_span = linspace(t0, t0+T, N_to_return);
    [t_sim,x_sim] = ode45(@two_body, t_span, x0);
    t = t_sim';
    x = x_sim';
end
end

function out = disturbance_model(t)
persistent generator
if t==0
    generator = RandStream('mt19937ar','Seed',0);
end
magnitude = 9.2e-6; % 1367 W/m2 * 2 m2 / 3e8 m/s
out = magnitude*(2*generator.rand(3,1) - 1);
if norm(out) > magnitude
    out = out*magnitude/norm(out);
end
end

function out = my_ode4(tspan, x0)
% This function performs runge-kutta integrations that always ends at the times where
% disturbances change.
persistent last_disturbance_time w
if tspan(1)==0
    last_disturbance_time = 0;
    w = disturbance_model(0);
end
step_disturbance = 3;
step_integration = 5; % ODE45 is using about a 10 second time step
out = zeros(length(x0), length(tspan));
t = tspan(1);
y = x0;
out(:,1) = x0;
next_req_index = 2;
while t < tspan(end)
    [t_next, select] = min([tspan(next_req_index), last_disturbance_time + step_disturbance, t + step_integration]);
    h = t_next - t;
    y = runge_kutta(t, y, h, w);
    if select==1
        out(:,next_req_index) = y;
        next_req_index = next_req_index + 1;
    end
    if t_next >= last_disturbance_time + step_disturbance
        w = disturbance_model(t_next);
        last_disturbance_time = t_next;
    end
    t = t_next;
end
end

function out = runge_kutta(t,y,h,u)
k1 = feval(@two_body,t,y,u);
k2 = feval(@two_body,t+h/2,y+k1*h/2,u);
k3 = feval(@two_body,t+h/2,y+k2*h/2,u);
k4 = feval(@two_body,t+h,y+h*k3,u);
out = y + h/6*(k1+2*k2+2*k3+k4);
end