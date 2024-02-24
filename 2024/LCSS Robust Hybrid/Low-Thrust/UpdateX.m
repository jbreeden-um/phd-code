function x1 = UpdateX(t0, x0, u, dt, use_disturbance)
if nargin==5 && use_disturbance
    t1 = t0 + dt;
    y = my_ode4([t0 t1], x0, u);
    x1 = y(:,2);
else
    t1 = t0 + dt;
    [~,x_sim] = ode45(@two_body, [t0 t1], x0, odeset('RelTol', 1e-10, 'AbsTol', 1e-10), u);
    x1 = x_sim(end,:)';
end
end

function out = disturbance_model_unforced(t)
global wc_max
persistent generator
if t==0
    generator = RandStream('mt19937ar','Seed',0);
end
out = generator.randn(3,1);
if norm(out) > 1
    out = out/norm(out);
end
out = wc_max*out;
end

function out = disturbance_model_forced(t, u_norm)
global wg_rate
persistent generator
if t==0
    generator = RandStream('mt19937ar','Seed',0);
end
direction = generator.randn(3,1);
direction = direction/norm(direction);
out = wg_rate*generator.rand*u_norm*direction;
end

function out = my_ode4(tspan, x0, u)
% This function performs runge-kutta integrations that always ends at the times where
% disturbances change.
persistent last_disturbance_time wc
if tspan(1)==0
    last_disturbance_time = 0;
    wc = disturbance_model_unforced(0);
end
step_disturbance = 300;
step_integration = 10;
% This time step leads to about 0.03 meters of error after 10 days, which I consider acceptable
% I am deliberately trading off speed of computation for accuracy, as the
% closed-loopedness should counteract any numerical innaccuracies.
out = zeros(length(x0), length(tspan));
t = tspan(1);
y = x0;
out(:,1) = x0;
next_req_index = 2;
wg = disturbance_model_forced(tspan(1), norm(u));
while t < tspan(end)
    [t_next, select] = min([tspan(next_req_index), last_disturbance_time + step_disturbance, t + step_integration]);
    h = t_next - t;
    y = runge_kutta(t, y, h, u + wg + wc);
    if select==1
        out(:,next_req_index) = y;
        next_req_index = next_req_index + 1;
    end
    if t_next >= last_disturbance_time + step_disturbance
        wc = disturbance_model_unforced(t_next);
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