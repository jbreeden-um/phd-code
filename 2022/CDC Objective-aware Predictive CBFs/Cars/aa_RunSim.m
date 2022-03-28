% Runs the simulation
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F); clear F;

global sim_case
sim_case = 2; % 1 or 2
control_case = 1; % 1, 2, or 3
verbose = 0; % recommended when sim_case=1, but this works in every case

dt = 0.01;
t = 0:dt:8;
N = length(t);
x0 = [-37; 10; -40; 10];

x = zeros(4,N)*NaN;
x(:,1) = x0;
u = zeros(2,N-1)*NaN;
H = zeros(1,N-1)*NaN;
h = zeros(1,N-1)*NaN;
compute = zeros(1,N-1)*NaN;

if verbose
    figure(8); clf;
    Hp = plot(t(1:end-1), H);
    xlabel 'Time (s)'; ylabel 'H^*';
    axis([0 7 -55 10]);

    figure(9); clf;
    up1 = plot(t(1:end-1), u(1,:)); hold on;
    up2 = plot(t(1:end-1), u(2,:));
    xlabel 'Time (s)'; ylabel 'u'; legend u_1 u_2;
    axis([0 7 -10 10]);

    figure(10); clf;
    hp = plot(t(1:end-1), h); hold on;
    plot([0 8], [0 0], 'r--');
    xlabel 'Time (s)'; ylabel 'Predicted Safety';
    axis([0 8 -60 10]);
end

for i=1:(N-1)
    h(i) = h_func(t(i),x(:,i));
    tic
    if control_case==1 || control_case==2
        [u(:,i), H(i)] = CalculateU(t(i),x(:,i),control_case);
    elseif control_case==3
        [u(:,i), H(i)] = CalculateU_MPC(t(i),x(:,i));
    else
        error('Unknown controller.');
    end
    compute(i) = toc;
    x(:,i+1) = UpdateX(t(i), t(i+1), x(:,i), u(:,i));
    
    if verbose
        t_pred = t(i):.1:(t(end)+1);
        h_pred = zeros(size(t_pred));
        for j=1:length(t_pred)
            h_pred(j) = h_func(t_pred(j), path_func(t_pred(j), t(i), x(:,i)));
        end
        set(hp, 'XData', t_pred, 'YData', h_pred);
        set(up1, 'YData', u(1,:));
        set(up2, 'YData', u(2,:));
        set(Hp, 'YData', H);
        drawnow;
    end
    
    waitbar(i/N);
end
mean_compute = mean(compute)

%%
figure(1); clf;
l1 = lane1(x0(1));
l2 = lane2(x0(3));
p1 = plot(l1(1), l1(2), 'bo', 'MarkerFaceColor', 'b'); hold on;
p2 = plot(l2(1), l2(2), 'go', 'MarkerFaceColor', 'g');
plot([-50, -3], [3, 3], 'k')
plot([-50, -3], [-3, -3], 'k')
plot([50, 3], [3, 3], 'k')
plot([50, 3], [-3, -3], 'k')
plot([-3, -3], [-50, -3], 'k')
plot([3, 3], [-50, -3], 'k')
plot([-3, -3], [50, 3], 'k')
plot([3, 3], [50, 3], 'k')
axis equal; axis([-50 50 -50 50]);
xlabel 'x (meters)';
ylabel 'y (meters)';
for i=1:N
    l1 = lane1(x(1,i));
    l2 = lane2(x(3,i));
    set(p1, 'XData', l1(1), 'YData', l1(2));
    set(p2, 'XData', l2(1), 'YData', l2(2));
    drawnow;
    if dt == 0.1
        pause(.05)
    end
end

figure(2); clf;
plot(t(1:end-1), u);
xlabel 'Time (s)';
ylabel 'u (m/s^2)';

figure(3); clf;
plot(t, x([1,3],:));
xlabel 'Time (s)';
ylabel 'z_i (meters)';
legend z_1 z_2

figure(4); clf;
plot(t, x([2,4],:)); hold on;
plot([t(1), t(end)], [12, 12], 'k--');
xlabel 'Time (s)';
ylabel '$$\dot{z}_i$$ (m/s)' interpreter latex;
legend({'$$\dot{z}_1$$','$$\dot{z}_2$$'},'interpreter','latex');

figure(5); clf;
plot(t(1:end-1), H); hold on;
plot(t(1:end-1), h);
xlabel 'Time (s)';
ylabel 'Constraints';
legend H^* h