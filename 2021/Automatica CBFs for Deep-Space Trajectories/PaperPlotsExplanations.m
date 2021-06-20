% Plots some of the expalantory figures
% Note that the switching figure is in the PaperPlotsFlyby3 scripts instead

%% Explanation of the restricted safe set
h = -1:.01:1;
hdot = -1:.01:1;
[x, y] = meshgrid(h, hdot);
a_max = 1;
absSq = @(x) x.*abs(x);
H = x + absSq(y)/(2*a_max);
figure(1); clf;
surf(x,y,H,'EdgeAlpha',0,'FaceColor','w'); 
hold on;
view([0;0;1]);

absRoot = @(x) sign(x)*sqrt(abs(x));
fxhatchfill(-1.5, 42, -3, 0, -1, 0, @(x) -1, @(x) min(1,absRoot(-2*a_max*x)), [0;1;0]);
fxhatchfill(1.5, 28, -1, 0.5, 0, 0.5, @(x) -1, @(x) min(1,absRoot(-2*a_max*x)), [1;0;0]);

h_crit = -absSq(hdot)/(2*a_max);
s_h = plot(h_crit,hdot,'k', 'LineWidth', 2); hold on;
s_com = plot([0 0], [-1 0], 'LineWidth', 3, 'Color', [0.8;0;0.8]);
axis equal;
axis([-1 1 -1 1 -2 2]);

xlabel('$h(t,x)$', 'interpreter', 'latex', 'FontSize', 13); 
ylabel('$\dot{h}(t,x)$', 'interpreter', 'latex', 'FontSize', 13);
legend([s_h, s_com],{'$\partial S_H$', '$S_H^\textrm{com}$'},'interpreter','latex','location','northeast','fontsize',13);
text(-0.6,-0.03,0.1,'$S_H^\textrm{res}$','interpreter','latex','fontsize',15);
text(0.03,-0.87,0.2,'$S_H^\textrm{bad}$','interpreter','latex','fontsize',15);
text(0.35,-0.03,1,'$\mathbf{R}^n \setminus S_H$','interpreter','latex','fontsize',15);

%% Explanation of the various maximizers
h0 = -1;
hdot0 = -3;
hddot0 = 2.2;
hdddot = -.4034;
t = 0:.1:5;
p = [hdddot, hddot0, hdot0, h0];
figure(2); clf;
plot(t, polyval(p,t*0.75),'b','LineWidth',2);
hold on;
p1 = plot(0,-1,'ro','LineWidth',2.5);
f = @(t) -polyval(p,t*0.75);
ts = fminsearch(f,3.6);
p2 = plot(ts,-1,'mo','LineWidth',2.5);
set(gcf,'Position',[2513 722 560 200]);
xlabel('Time $\beta$ since $t$', 'interpreter','latex','fontsize',13);
ylabel('$h(t+\beta,\chi(\beta,t,x))$','interpreter','latex','fontsize',13);

plot(0.8*[1 1],[-3 -1],'k--');
plot(2.9*[1 1],[-3 -1],'k--');
plot(4.25*[1 1],[-3 -1],'k--');
text(0.08,-2.75,'$\mathcal{B}_{loc}(0)$', 'interpreter','latex','fontsize',15)
text(3.15,-2.75,'$\mathcal{B}_{loc}(3.6)$', 'interpreter','latex','fontsize',15)

function fxhatchfill(slope, n, xstart, xend, xmin, xmax, ybot, ytop, color)
for i=linspace(xstart, xend, n)
    x0 = i;
    if slope > 0
        ys = ybot;
        yf = ytop;
    elseif slope < 0
        ys = ytop;
        if x0 > -0.5
            yf = ytop;
        else
            yf = ybot;
        end
    else
        error('Horizontal hatch not allowed');
    end
    y0 = ys(x0);
    if x0 < xmin
        y_left = ys(x0) + slope*(xmin - x0);
        if y_left > ybot(xmin) && y_left < ytop(xmin)
            x0 = xmin;
            y0 = y_left;
        else
            continue;
        end
    end
    f = @(x) y0 + slope*(x - x0) - yf(x);
    x1 = fzero(f, xend);
    y1 = yf(x1);
    if x1 > xmax
        x1 = xmax;
        y1 = y0 + slope*(xmax - x0);
    end
    plot([x0, x1], [y0, y1], 'Color', color);
end
end
