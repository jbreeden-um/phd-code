% Runs the simulation
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F); clear F;

global debris_ke time_step
control_case = 0; % 0, 1, or 2
verbose = 0; % recommended when control_case=1, but this works in every case

mu = 398600;
dt = 5;
t = 0:dt:2500;
N = length(t);
x0 = KeplerToCartesian(7000, 0, deg2rad(23.4), pi/2, 0, deg2rad(10)+1, mu);
debris_ke = struct('a', 7000, 'e', 0, 'i', deg2rad(113.4), 'O', 3*pi/2, 'o', 0, 'nu0', deg2rad(190)+1);

time_step = 0.01;
t_all = 0:time_step:t(end);
N_per_step = dt/time_step;

x = zeros(6,(N-1)*N_per_step+1)*NaN;
x(:,1) = x0;
u = zeros(3,N-1)*NaN;
H = zeros(1,N-1)*NaN;
h = zeros(1,N-1)*NaN;
compute = zeros(1,N-1)*NaN;

if verbose
    figure(8); clf;
    Hp = plot(t(1:end-1), H);
    xlabel 'Time (s)'; ylabel 'H^*';

    figure(9); clf;
    up1 = plot(t(1:end-1), u(1,:)); hold on;
    up2 = plot(t(1:end-1), u(2,:));
    up3 = plot(t(1:end-1), u(3,:));
    xlabel 'Time (s)'; ylabel 'u'; legend u_1 u_2;

    figure(10); clf;
    hp = plot(t(1:end-1), H); hold on;
    plot([0 t(end)], [0 0], 'r--');
    xlabel 'Time (s)'; ylabel 'Predicted Safety';
end

for i=1:(N-1)
    x_curr = x(:,1+N_per_step*(i-1));
    h(i) = h_func(t(i), x_curr);
    tic
    if control_case==1 || control_case==2
        [u(:,i), H(i)] = CalculateU(t(i),x_curr,control_case);
    elseif control_case==0
        u(:,i) = [0;0;0];
    else
        error('Unknown controller.');
    end
    compute(i) = toc;
    x(:,(1+N_per_step*(i-1)):(1+N_per_step*i)) = UpdateX(t(i), t(i+1), x_curr, u(:,i))';
    
    if verbose && mod(t(i),25) == 0
        t_pred = t(i):10:(t(end)+1);
        h_pred = zeros(size(t_pred));
        for j=1:length(t_pred)
            h_pred(j) = h_func(t_pred(j), path_func(t_pred(j), t(i), x_curr));
        end
        set(hp, 'XData', t_pred, 'YData', h_pred);
        set(up1, 'YData', u(1,:));
        set(up2, 'YData', u(2,:));
        set(up3, 'YData', u(3,:));
        set(Hp, 'YData', H);
        drawnow;
    end
    
    waitbar(i/N);
end
mean_compute = mean(compute)

%%
check = h > max(h)-100;
[~, index1] = max(check);
[~, index2] = min(check(index1:end));
index2 = index2 + index1 - 1;
t_extra = t(index1):time_step:t(index2-1);
t_dist = [t(1:index1-1), t_extra, t(index2:end-1)];
h_dist = [h(1:index1-1), NaN*t_extra, h(index2:end)];
[~, offset] = min(abs(t_all - t_dist(index1)));
for i=index1:(index1+length(t_extra)-1)
    index = offset+(i-index1);
    if abs(t_all(index) - t_dist(i)) > time_step^2
        error('Unmmatched indexes');
    end
    h_dist(i) = h_func(t_dist(i), x(:,index));
    waitbar((i - index1)/(index1+length(t_extra)-1-index1));
end

%%
figure(2); clf;
plot(t(1:end-1), u);
xlabel 'Time (s)';
ylabel 'u (km/s^2)';

figure(3); clf;
plot(t_all, x(1:3,:));
xlabel 'Time (s)';
ylabel 'r (km)';
legend r_1 r_2 r_3

figure(4); clf;
plot(t_all, x(4:6,:)); hold on;
plot([t(1), t(end)], [12, 12], 'k--');
xlabel 'Time (s)';
ylabel 'v (km/s)';
legend v_1 v_2 v_3

figure(5); clf;
plot(t(1:end-1), H); hold on;
plot(t_dist, h_dist);
xlabel 'Time (s)';
ylabel 'Constraints';
legend H^* h