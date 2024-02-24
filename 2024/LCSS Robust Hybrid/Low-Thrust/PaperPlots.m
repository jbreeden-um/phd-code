load sim_data;

figure(1); clf;
scale = 1e3;
[~, vertices, indices] = make_die;
r = norm(mean(vertices(:,indices(:,1)),2));
rho = 1e4;
vertices = vertices*rho/r;
trimesh(indices',vertices(1,:)/scale,vertices(2,:)/scale,vertices(3,:)/scale,'EdgeColor','k','FaceColor',[0.4;0.4;0.4],'FaceAlpha',0); hold on;
xc = zeros(size(x));
for i=1:length(x)
    xc(:,i) = get_center(t(i));
end
dx = x - xc;
dxhat = xhat - xc;
p1 = plot3(dx(1,:)/scale, dx(2,:)/scale, dx(3,:)/scale, 'b', 'LineWidth', 1);
p2 = plot3(dxhat(1,:)/scale, dxhat(2,:)/scale, dxhat(3,:)/scale, 'Color', [0; 0.7; 0], 'LineWidth', 1);
axis equal;
xlabel 'x (km)'; ylabel 'y (km)'; zlabel 'z (km)';
legend([p1 p2],{'Real','Estimate'});
set(gcf, 'Position', [50 550 560 420]);
view(-200,6)

try
    hhat_center(1);
catch
    hhat_center = zeros(size(hhat));
    for i=1:i_end
        [~,~,~,hhat_center(:,i)] = CBF_icosa(t(i),xhat(:,i),rhohat(:,i),u(:,i));
        waitbar(i/i_end);
    end
end
figure(2); clf;
td = t/(24*3600);
tt = [td, fliplr(td)];
yy = [hhat_max + hhat_error, fliplr(hhat_max - hhat_error)];
jj = ~isnan(yy);
fill(tt(jj), yy(jj)/scale, 'r', 'FaceAlpha', 0.6, 'EdgeAlpha', 0); hold on;
h_max = max(h,[],1);
hhat_max = max(hhat_center,[],1);
hhat_error = max(hhat,[],1) - hhat_max;
plot(td, h_max/scale, 'b'); hold on;
% plot(td, hhat_max/scale, 'r');
set(gcf, 'Position', [50 250 600 200]);
axis([0 5 -12 0]);
xlabel 'Time (days)'; ylabel 'max(h_i)  (km)';