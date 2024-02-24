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
Kdd_f = @(x1, xc) (2*(xc(1) - x1(1))*(xc(3) - x1(3)) + 2*(xc(2) - x1(2))*(xc(4) - x1(4)))^2/(4*((xc(1) - x1(1))^2 + (xc(2) - x1(2))^2)^(3/2)) + (2*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2))*(xc(1) - x1(1)) + 2*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2))*(xc(2) - x1(2)) - 2*(xc(3) - x1(3))^2 - 2*(xc(4) - x1(4))^2)/(2*((xc(1) - x1(1))^2 + (xc(2) - x1(2))^2)^(1/2));
period = 2*pi*sqrt(oe.a^3/mu);

% The following x0 is used to compute Kdd_max_arr
center_func = @get_center;

% The following x0 is used to compute Kdd_compare. The values of Kdd_max for different
% obstacle indices are slightly different, but this difference is very slight.
% for i=1:6, obstacle_location(i,0); end
% center_func = @(t) obstacle_location(6, t);

level = 5e7;
rel_pos = sqrt(level/min(eig(P)));
rel_vel = sqrt(level/max(eig(P)));
b = [rel_pos; rel_pos; rel_vel; rel_vel]*1.2;

obj = @(x) -Kdd_f(x(1:4), center_func(x(5)));
N = 100;
x_crit = [0;0;0;0;0];
f_crit = Inf;
f_crit_all = [0,0,0];
x_crit_all = zeros(5,3);
tic
for j=0:2
    t_guess = period*j/2;
    for i=1:N
        x_guess = [center_func(t_guess) + randn(4,1).*b/1.5; t_guess];
        [x_val, f_val] = fmincon(obj, x_guess, [], [], [], [], [-1e7;-1e7;-1e7;-1e7; -0.05*period], [1e7;1e7;1e7;1e7; 1.05*period], ...
            @(x) nonlin(x(1:4), P, level, center_func(x(5))), optimset('Display', 'off'));
%         [x_val, f_val] = fmincon(@(x) obj([x;t_guess]), x_guess(1:4), [], [], [], [], center_func(t_guess)-b, center_func(t_guess)+b, ...
%             @(x) nonlin(x(1:4), P, level, center_func(t_guess)), optimset('Display', 'off'));
        if f_val < f_crit
            f_crit = f_val;
            x_crit = x_val;
        end
%         waitbar((j*N+i)/(3*N));
        waitbar(i/N);
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

% Sample values of f_crit_all for obstacle 6
% f_crit_all =
%   -0.005355287955695  -0.002886495133254  -0.005351811083131
% In conclusion, focusing on apoapsis or periapsis provides the correct number. 
% The results are close to the circular results, so this all makes sense.

%%
levels = [5e4 5e5 5e6 5e7];
levels_sq = sqrt(levels);
Kdd_max_arr = [1.641061759398710e-04 5.253583594431177e-04 0.001684155744250 0.005353582703918];
plot(levels_sq, Kdd_max_arr);

Kdd_compare = [0.005352080239788; 0.005350817362626; 0.005355339262621; 0.005322748674985; 0.005250501552609; 0.005355287955695]; 
compare = Kdd_compare / Kdd_max_arr(4)

m = max(Kdd_max_arr./levels_sq)*max(compare);
% hddot_max vs sqrt(V) follows an almost exactly linear curve with the above coefficient
% Encode m as 7.6e-7 to be conservative

function [C, Ceq] = nonlin(x, P, level, x0)
x1 = x(1:4);
C = (x1 - x0)'*P*(x1 - x0) - level;
Ceq = [];
end