function PaperPlots
%%
sim1 = load('data/sim_data120');
sim2 = load('data/sim_data240');
sim3 = load('data/sim_data360');

figure(1); clf;
theta = linspace(0, 2*pi, 100);
cx = cos(theta);
cy = sin(theta);
rho = 1e3;
scale = 1e3;
t_obs = 100;
obs = convert_to_vis(t_obs, obstacle_location(1,t_obs)); fill((obs(1)+cx*rho)/scale, (obs(2)+cy*rho)/scale, 'r'); hold on;
obs = convert_to_vis(t_obs, obstacle_location(2,t_obs)); fill((obs(1)+cx*rho)/scale, (obs(2)+cy*rho)/scale, 'r');
obs = convert_to_vis(t_obs, obstacle_location(3,t_obs)); fill((obs(1)+cx*rho)/scale, (obs(2)+cy*rho)/scale, 'r');
obs = convert_to_vis(t_obs, obstacle_location(4,t_obs)); fill((obs(1)+cx*rho)/scale, (obs(2)+cy*rho)/scale, 'r');
obs = convert_to_vis(t_obs, obstacle_location(5,t_obs)); fill((obs(1)+cx*rho)/scale, (obs(2)+cy*rho)/scale, 'r');
obs = convert_to_vis(t_obs, obstacle_location(6,t_obs)); fill((obs(1)+cx*rho)/scale, (obs(2)+cy*rho)/scale, 'r');
fill([-30 30 30 -30 -30], [0 0 10 10 0], 'r');

xhat1a = sim1.xhat_lin(1:2,:) + sim1.rhohat(1,:);
xhat1b = sim1.xhat_lin(1:2,:) - sim1.rhohat(1,:);
fill([xhat1a(1,:), fliplr(xhat1b(1,:))]/scale, [xhat1a(2,:), fliplr(xhat1b(2,:))]/scale, 'b', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.1)
xhat2a = sim2.xhat_lin(1:2,:) + sim2.rhohat(1,:);
xhat2b = sim2.xhat_lin(1:2,:) - sim2.rhohat(1,:);
fill([xhat2a(1,:), fliplr(xhat2b(1,:))]/scale, [xhat2a(2,:), fliplr(xhat2b(2,:))]/scale, 'b', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.1)
xhat3a = sim3.xhat_lin(1:2,:) + sim3.rhohat(1,:);
xhat3b = sim3.xhat_lin(1:2,:) - sim3.rhohat(1,:);
fill([xhat3a(1,:), fliplr(xhat3b(1,:))]/scale, [xhat3a(2,:), fliplr(xhat3b(2,:))]/scale, 'b', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.1)

p1 = plot(sim1.x_lin(1,:)/scale, sim1.x_lin(2,:)/scale, 'b', 'LineWidth', 2);
p2 = plot(sim2.x_lin(1,:)/scale, sim2.x_lin(2,:)/scale, 'b--', 'LineWidth', 2);
p3 = plot(sim3.x_lin(1,:)/scale, sim3.x_lin(2,:)/scale, 'b:', 'LineWidth', 2);

legend([p1 p2 p3], {'120 s', '240 s', '360 s'}, 'FontSize', 20, 'Location', 'SouthEast');

axis equal;
axis([-20 20 -31 2]);
set(gcf, 'Position', [50 550 560 420]);
set(gca, 'FontSize', 15, 'XTick', -20:5:20);
xlabel('x_1 (km)','FontSize',18); ylabel('x_2 (km)','FontSize',18);

%%
hhat_center = zeros(7, length(sim3.t))*NaN;
for i=1:sim3.i_end
    [~,~,~,hhat_center(1,i)] = CBF_obs(sim3.t(i),sim3.xhat(:,i),sim3.TM,[0;0],[],1,'post');
    [~,~,~,hhat_center(2,i)] = CBF_obs(sim3.t(i),sim3.xhat(:,i),sim3.TM,[0;0],[],2,'post');
    [~,~,~,hhat_center(3,i)] = CBF_obs(sim3.t(i),sim3.xhat(:,i),sim3.TM,[0;0],[],3,'post');
    [~,~,~,hhat_center(4,i)] = CBF_obs(sim3.t(i),sim3.xhat(:,i),sim3.TM,[0;0],[],4,'post');
    [~,~,~,hhat_center(5,i)] = CBF_obs(sim3.t(i),sim3.xhat(:,i),sim3.TM,[0;0],[],5,'post');
    [~,~,~,hhat_center(6,i)] = CBF_obs(sim3.t(i),sim3.xhat(:,i),sim3.TM,[0;0],[],6,'post');
    [~,~,~,hhat_center(7,i)] = CBF_dock(sim3.t(i),sim3.xhat(:,i),sim3.TM,[0;0],[],'post');
    waitbar(i/sim3.i_end);
end

hhat_error = sim3.hhat - hhat_center;
colors = [                   0   0.447000000000000   0.741000000000000;
          0.850000000000000   0.325000000000000   0.098000000000000;
          0.929000000000000   0.694000000000000   0.125000000000000;
          0.494000000000000   0.184000000000000   0.556000000000000;
          0.466000000000000   0.674000000000000   0.188000000000000;
          0.301000000000000   0.745000000000000   0.933000000000000;
          0.635000000000000   0.078000000000000   0.184000000000000];

figure(2); clf;
for i=1:7
    fill([sim3.t, fliplr(sim3.t)], [hhat_center(i,:) + hhat_error(i,:), fliplr(hhat_center(i,:) - hhat_error(i,:))]/scale, ...
        colors(i,:), 'FaceAlpha', 0.4, 'EdgeAlpha', 0.1);
    hold on;
end
for i=1:7
    plot(sim3.t, sim3.h(i,:)/scale, 'Color', colors(i,:));
end
axis([0 5000 -35 0]);
xlabel 'Time (seconds)'; ylabel 'h_i (km)';
set(gcf, 'Position', [50 250 600 165]);

end