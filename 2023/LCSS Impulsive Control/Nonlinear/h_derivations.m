% This file computes \ddot{\kappa}_{max} for use in the CBF_obs file

syms mu real
syms x11(t) x12(t) x13(t) x14(t) x01(t) x02(t) x03(t) x04(t)
x1d = [x13; x14; -mu*[x11; x12]/(x11^2 + x12^2)^(3/2)];
x0d = [x03; x04; -mu*[x01; x02]/(x01^2 + x02^2)^(3/2)];

K = -sqrt((x11-x01)^2 + (x12-x02)^2);
Kd = subs(diff(K), [diff(x11); diff(x12); diff(x13); diff(x14); diff(x01); diff(x02); diff(x03); diff(x04)], [x1d; x0d])
Kdd = subs(diff(Kd), [diff(x11); diff(x12); diff(x13); diff(x14); diff(x01); diff(x02); diff(x03); diff(x04)], [x1d; x0d])

%%
exp = string(Kdd);
exp = replace(exp, {'x11(t)', 'x12(t)', 'x13(t)', 'x14(t)'}, {'x1(1)', 'x1(2)', 'x1(3)', 'x1(4)'});
exp = replace(exp, {'x01(t)', 'x02(t)', 'x03(t)', 'x04(t)'}, {'xc(1)', 'xc(2)', 'xc(3)', 'xc(4)'});
fprintf('\n%s\n\n',exp);

%%
mu = 398600e9;
Re = 6378e3;
alt = 600e3;
n = sqrt(mu/(Re+alt)^3);
A_sys = [0,     0, 1,    0;
         0,     0, 0,    1;
         3*n^2, 0, 0,    2*n;
         0,     0, -2*n, 0];
B_sys = [zeros(2); eye(2)];
R = eye(2)*100;
Q = eye(4);
[~, P, ~] = lqr(A_sys,B_sys,Q,R);
Kdd_f = @(x1, xc) (2*(xc(1) - x1(1))*(xc(3) - x1(3)) + 2*(xc(2) - x1(2))*(xc(4) - x1(4)))^2/(4*((xc(1) - x1(1))^2 + (xc(2) - x1(2))^2)^(3/2)) + (2*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2))*(xc(1) - x1(1)) + 2*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2))*(xc(2) - x1(2)) - 2*(xc(3) - x1(3))^2 - 2*(xc(4) - x1(4))^2)/(2*((xc(1) - x1(1))^2 + (xc(2) - x1(2))^2)^(1/2));

% The following x0 is used to compute Kdd_max_arr
x0 = get_center(0);
% The following x0 is used to compute Kdd_compare. The resultant values of Kdd_max are
% slightly different depending on the obstacle used, but this difference is very slight.
% [r0, v0] = obstacle_location(4,0); x0 = [r0; v0];

level = 5e7;
rel_pos = sqrt(level/min(eig(P)));
rel_vel = sqrt(level/max(eig(P)));
b = [rel_pos; rel_pos; rel_vel; rel_vel]*1.1;

obj = @(x) -Kdd_f(x, x0);
N = 100;
x_crit = [0;0;0;0];
f_crit = Inf;
for i=1:N
    x_guess = x0 + randn(4,1).*b/1.5;
    [x_val, f_val] = fmincon(obj, x_guess, [], [], [], [], x0-b, x0+b, @(x) nonlin(x, P, level, x0), optimset('Display', 'off'));
    if f_val < f_crit
        f_crit = f_val;
        x_crit = x_val;
    end
    waitbar(i/N);
end
Kdd_max = -f_crit
check_negative = nonlin(x_crit, P, level, x0)

%%
levels = [5e4 5e5 5e6 5e7];
levels_sq = sqrt(levels);
Kdd_max_arr = [1.597610127292981e-04 5.132459527139588e-04 0.001680074902524 0.005348439517149];
plot(levels, Kdd_max_arr);

Kdd_compare = [0.005350735588394; 0.005349076301030; 0.005186462209707; 0.005227767886254]; 
compare = Kdd_compare / Kdd_max_arr(4)

m = 7.586e-07*max(compare);
% hddot_max vs sqrt(V) follows an almost exactly linear curve with the above coefficient

function [C, Ceq] = nonlin(x, P, level, x0)
x1 = x(1:4);
C = (x1 - x0)'*P*(x1 - x0) - level;
Ceq = [];
end