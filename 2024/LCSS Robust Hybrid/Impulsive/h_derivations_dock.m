% This file computes \ddot{\kappa}_{max} for use in the CBF_dock file. Because the docking
% CBF is a planar exclusion rather than a circular exclusion, the value of
% \ddot{\kappa}_{max} is much larger.

syms mu real
syms x11(t) x12(t) x13(t) x14(t) x01(t) x02(t) x03(t) x04(t) v
x1d = [x13; x14; -mu*[x11; x12]/(x11^2 + x12^2)^(3/2)];
x0d = [x03; x04; -mu*[x01; x02]/(x01^2 + x02^2)^(3/2)];

% v is the velocity of the target
% For a circular orbit, v is constant. For an elliptical one, v is time-varying
% However, for simplicity, I am going to let v be a scaling constant rather than the 
% instantaneous value of the real velocity.
K = ( (x11-x01)*x03 + (x12-x02)*x04 ) / v;
Kd = subs(diff(K), [diff(x11); diff(x12); diff(x13); diff(x14); diff(x01); diff(x02); diff(x03); diff(x04)], [x1d; x0d])
Kdd = subs(diff(Kd), [diff(x11); diff(x12); diff(x13); diff(x14); diff(x01); diff(x02); diff(x03); diff(x04)], [x1d; x0d])

%%
exp = string(Kdd);
exp = replace(exp, {'x11(t)', 'x12(t)', 'x13(t)', 'x14(t)'}, {'x1(1)', 'x1(2)', 'x1(3)', 'x1(4)'});
exp = replace(exp, {'x01(t)', 'x02(t)', 'x03(t)', 'x04(t)'}, {'xc(1)', 'xc(2)', 'xc(3)', 'xc(4)'});
fprintf('\n%s\n\n',exp);

%%
mu = 398600e9;
[~, ~, oe] = get_center(0);
n = sqrt(mu/oe.a^3);
% This linear system does not matter at all. It is just a way to generate a P matrix.
A_sys = [0,     0, 1,    0;
         0,     0, 0,    1;
         3*n^2, 0, 0,    2*n;
         0,     0, -2*n, 0];
B_sys = [zeros(2); eye(2)];
R = eye(2)*100;
Q = eye(4);
[~, P, ~] = lqr(A_sys,B_sys,Q,R);
v = sqrt(mu/oe.a); % scaling constant, not important
Kdd_f = @(x1, xc) (xc(3)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + xc(4)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2)) + (mu*xc(3)*(xc(1) - x1(1)))/(xc(1)^2 + xc(2)^2)^(3/2) + (2*mu*xc(1)*(xc(3) - x1(3)))/(xc(1)^2 + xc(2)^2)^(3/2) + (mu*xc(4)*(xc(2) - x1(2)))/(xc(1)^2 + xc(2)^2)^(3/2) + (2*mu*xc(2)*(xc(4) - x1(4)))/(xc(1)^2 + xc(2)^2)^(3/2) - (3*mu*xc(1)*(xc(1) - x1(1))*(2*xc(1)*xc(3) + 2*xc(2)*xc(4)))/(2*(xc(1)^2 + xc(2)^2)^(5/2)) - (3*mu*xc(2)*(xc(2) - x1(2))*(2*xc(1)*xc(3) + 2*xc(2)*xc(4)))/(2*(xc(1)^2 + xc(2)^2)^(5/2)))/v;

period = 2*pi*sqrt(oe.a^3/mu);

level = 5e4;
rel_pos = sqrt(level/min(eig(P)));
rel_vel = sqrt(level/max(eig(P)));
b = [rel_pos; rel_pos; rel_vel; rel_vel]*1.1;

obj = @(x) -Kdd_f(x(1:4), get_center(x(5)));
N = 100;
x_crit = [0;0;0;0;0];
f_crit = Inf;
f_crit_all = [0,0,0];
x_crit_all = zeros(5,3);
tic
for j=0:2
    t_guess = period*j/2;
    if j==1, t_guess = rand*period; end
    for i=1:N
        x_guess = [get_center(t_guess) + randn(4,1).*b/1.5; t_guess];
        [x_val, f_val] = fmincon(obj, x_guess, [], [], [], [], [-1e7;-1e7;-1e7;-1e7; -0.05*period], [1e7;1e7;1e7;1e7; 1.05*period], ...
            @(x) nonlin(x(1:4), P, level, get_center(x(5))), optimset('Display', 'off'));
        if f_val < f_crit
            f_crit = f_val;
            x_crit = x_val;
        end
        waitbar((j*N+i)/(3*N));
    end
    f_crit_all(j+1) = f_crit;
    if length(x_crit)==5
        x_crit_all(:,j+1) = x_crit;
    else
        x_crit_all(:,j+1) = [x_crit; t_guess];
    end
    f_crit = Inf;
end
toc
Kdd_max = -min(f_crit_all)
if length(x_crit)==5
    check_negative = nonlin(x_crit(1:4), P, level, get_center(x_crit(5)))
else
    check_negative = nonlin(x_crit(1:4), P, level, get_center(t_guess))
end

%%
levels = [5e4 5e5 5e6 5e7];
levels_sq = sqrt(levels);
Kdd_max_arr = [0.104361153591622 0.330019076974629 1.043611662401607 3.300190579949684];
plot(levels_sq, Kdd_max_arr);

m = 0.000467;
% hddot_max vs sqrt(V) follows an almost exactly linear curve with the above coefficient (rounded up)

function [C, Ceq] = nonlin(x, P, level, x0)
x1 = x(1:4);
C = zeros(2,1);
C(1) = (x1 - x0)'*P*(x1 - x0) - level;
C(2) = dot(x1(1:2)-x0(1:2), x0(3:4)/norm(x0(3:4)));
Ceq = [];
end