%% Mass Parameter Setup
if 0
scales = 0:.4:3;
M2 = zeros(size(scales));
M3 = zeros(size(scales));
for i=1:length(scales)
    [M2(i), M3(i)] = get_M(scales(i));
    i
end
save('more_data1', 'scales', 'M2', 'M3');
end
%%
if 0
M3 = 2.738e-05*scales +  0.006144;
% dt = (.25:.625:5.5)*.2;
dt = linspace(0.05,1,8);
mu = zeros(length(dt), length(scales));
for it=1:length(dt)
    for im=1:length(scales)
        mu(it, im) = get_mu(dt(it), M2(im), M3(im));
        p = ((it-1)*length(scales)+im);
        progress = p/(length(dt)*length(scales))
        save(['Results/Temp/more_data2_' num2str(p)], 'scales', 'M2', 'M3', 'dt', 'mu');
    end
end
save('more_data2', 'scales', 'M2', 'M3', 'dt', 'mu');
end

%%
constants = load('parameters');
torque_axis = scales*constants.torque;
surf(torque_axis, dt, mu)
xlabel 'torque'; ylabel 'dt';

[xData, yData, zData] = prepareSurfaceData( torque_axis, dt, mu );
ft = fittype( 'b*x + a*x*y + d*y + e + g*y^2', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.346460985400705 0.37587321694393 0.225067882127913 0.350013919000766 0.287084613146176];
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

N = 200;
rx = linspace(torque_axis(1), torque_axis(end), N);
ry = linspace(dt(1), 0.32, 1.2*N);
M3z = 2.738e-05*rx/constants.torque +  0.006144;
M2z = constants.Z(3,3)*rx;
z = fitresult(rx', ry);
z0 = z;
z0(z < 0) = 0;
zcond = 2*M2z + M3z.*ry';
cutoff = min(min(z(z >= zcond')));
z0(z < zcond') = cutoff;

%%
figure(2); clf;
surf(ry, rx*1e6, z0, 'EdgeAlpha', 0);
ylabel '\xi_{max} (\muN-m)';
xlabel 'T (s)';
zlabel '\mu'
set(gca, 'FontSize', 14);

zcond = min(1.1*z', zcond);
m = max(z(:));
hold on;
% contour(ry, rx*1e6, z+m, [m m], 'LineWidth', 3, 'Color', 'k');
% surf(ry, rx*1e6, zcond', 'FaceColor', 'k', 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
surf(ry, rx*1e6, ones(size(z))*cutoff, 'FaceColor', 'k');

axis tight;
h = colorbar;
% ylabel(h, '\mu^*', 'FontSize', 14);
set(h, 'Ticks', 0:.0001:.0024)
view([0 0 1]);
% view(180,90)

plot3(0.2, constants.torque*1e6, m, 'wx', 'MarkerSize', 10, 'LineWidth', 3);

set(gcf, 'Position', [1200 600 560 350]);
set(gca, 'XTick', 0.05:0.05:32, 'XLim', [0.042 0.32], 'YLim', [-1, rx(end)*1e6]);
grid off;

function [M2plus, M3plus] = get_M(torque_scale)
%%
global ctheta
ctheta = cos(pi/4);
    % results are the same for any angle less than pi/2, and are overconservative for
    % angles larger than pi/2
tetrahedron_points = [0.5, 0, 0;
                      -0.5, 0, 0;
                      0, sqrt(3)/2, 0;
                      0, 0, 0];
f = @(x) [norm(x - tetrahedron_points(1,:)')^2 - 1;
          norm(x - tetrahedron_points(2,:)')^2 - 1;
          norm(x - tetrahedron_points(3,:)')^2 - 1];
x = fsolve(f, [0; 0.5; 0.5], optimset('Display', 'off'));
tetrahedron_points(4,:) = x';
global wheel_axis
wheel_axis(:,1) = cross(tetrahedron_points(1,:) - tetrahedron_points(2,:), tetrahedron_points(1,:) - tetrahedron_points(3,:))';
wheel_axis(:,2) = cross(tetrahedron_points(1,:) - tetrahedron_points(4,:), tetrahedron_points(1,:) - tetrahedron_points(2,:))';
wheel_axis(:,3) = cross(tetrahedron_points(1,:) - tetrahedron_points(3,:), tetrahedron_points(1,:) - tetrahedron_points(4,:))';
wheel_axis(:,4) = cross(tetrahedron_points(3,:) - tetrahedron_points(4,:), tetrahedron_points(2,:) - tetrahedron_points(3,:))';
    % the above has been ordered so that all the wheels positive spin axes face away from 
    % the center of the spacecraft
wheel_axis = wheel_axis./vecnorm(wheel_axis);
wheel_limit = 0.001; % https://www.cubespace.co.za/products/adcs-components/cubewheel/#cubewheel-specifications
wheel_limit = 0.7*wheel_limit; % to compensate for there being four wheels
wheel_momentum = 0.01082;
wheel_rate_limit = 2*pi*6000/60; % the max wheel rate in radians per second, note this leads to poor numerical conditioning
wheel_J = wheel_momentum / wheel_rate_limit;

global Z Jtot Jw
mass = 12; % kg, 6U CubeSat
d = [0.1; 0.2263; 0.3405]; % https://www.isispace.nl/product/6-unit-cubesat-structure/
Jb = 1/12*mass*diag([d(2)^2+d(3)^2; d(1)^2+d(3)^2; d(1)^2+d(2)^2]);
    % Assume a principal axis frame, since the wheels aren't a simple cube anyways.
Jw = diag([1;1;1;1])*wheel_J;
Jtot = Jb;
for i=1:4, Jtot = Jtot + wheel_J*wheel_axis(:,i)*wheel_axis(:,i)'; end
Z = inv([Jtot, wheel_axis*Jw; Jw*wheel_axis', Jw]);  
Je = inv(Z(1:3,1:3));

w_max = pi/180; % 1 deg/s absolute limit on the largest axis
E_max = [1,0,0]*Je*[1;0;0]*w_max^2;

A = d(2)*d(3)*2;
rho = 5.5e-12; % From the model in 42 (for fairly high solar activity)
v = sqrt(398600/(6378+500))*1e3; % 500 km altitude
Cd = 2.4; % guess from personal experience
lever = 1/2*d(3);
force = 1/2*Cd*A*rho*v^2;
torque = force*lever*torque_scale;

global O_Earth_from_Sol omega_Earth
O_Earth_from_Sol = [1, 0, 0; 0, cosd(-23.4392911), sind(-23.4392911); 0, -sind(-23.4392911), cosd(-23.4392911)]; % Ref: https://github.com/ericstoneking/42/blob/abd2b0bfbee8c10bf86e07d726340ef0bf0f1308/Source/42init.c
    % Orientation matrix converts vectors from ecliptic frame to Earth polar frame
omega_Earth = 2*pi/(365.25*24*3600);

year = 366*24*3600; % this is just used for scaling the time variable
N_max = 500;

%% Find M_2
% tic
% g1 = [zeros(3,4); zeros(3,4); eye(4)];
% g2 = [zeros(3,3); Z(1:3,1:3); zeros(4,3)];

% x(1) = t, x(2:4) = p, x(5:7) = omega, x(8:10) = d
% func = @(x) rgrad_hdot(x(1)*year, x(2:7))*g2*x(8:10);
p_con = @(x) rh(x(1)*year, x(2:7)); % in the safe set of quaternions
w_con = @(x) (x(5:7)'*Je*x(5:7) - E_max)*1e3; % in the safe set of angular velocities
d_con = @(x) norm(x(8:10)) - torque; % disturbance within limit
r_con = @(x) norm(x(2:4)) - 1; 
    % This should be an equality constraint, but the code is more efficient as is.
    % Note that later on we will enforce this equality, but here the worst case occurs
    % with equality anyways.
% lower = -[0; ones(3,1); 3*w_max*ones(3,1); torque*ones(3,1)];
% upper =  [1; ones(3,1); 3*w_max*ones(3,1); torque*ones(3,1)];
% nonlcon = @(x) assemble_constraint(x, @(y) [p_con(y); w_con(y); d_con(y); r_con(y)]);

M2minus = -Z(3,3)*torque;
% toc

%%
M2plus = Z(3,3)*torque;

%% Find M_3
% tic
% x(1) = t, x(2:4) = p, x(5:7) = omega, x(8:10) = d, x(11:14) = W, x(15:18) = wheel_commands
rng(0);
func = @(x) rphidot(x(1)*year, [x(2:7); x(11:14)*1e3], x(15:18), x(8:10));
lower = -[0; ones(3,1); w_max*ones(3,1)*3; torque*ones(3,1); wheel_rate_limit*ones(4,1)/1e3; wheel_limit*ones(4,1)];
upper =  [1; ones(3,1); w_max*ones(3,1)*3; torque*ones(3,1); wheel_rate_limit*ones(4,1)/1e3; wheel_limit*ones(4,1)];
nonlcon = @(x) assemble_constraint(x, @(y) [p_con(y); w_con(y); d_con(y); r_con(y)]);
N = 100;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(18,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon, optimset('Display','off'));
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
% min_val
if ~isempty(min_res)
    nonlcon(min_res);
    M3minus = min_val;
    M3plus = -min_val;
else
   M3plus = NaN;
   M3minus = NaN;
   disp('No solution to M3');
end
% toc
end

function res = get_mu(dt, M2, M3)
mu_test = .0002:.0002:.0024;
results = zeros(size(mu_test))*NaN;
count = 0;
min_guess = 1;
max_guess = length(results);
while count < length(mu_test) && sum(isnan(results)) > 0
    guess = ceil((min_guess+max_guess)/2);
%     guess = floor((min_guess+max_guess)/2);
    if isnan(results(guess))
        results(guess) = test_mu(mu_test(guess), dt, M2, M3);
    elseif isnan(results(guess-1))
        guess = guess-1;
        results(guess) = test_mu(mu_test(guess), dt, M2, M3);
    else
        return;
    end
    
    if results(guess) == 0
        res = mu_test(guess);
        return;
    elseif results(guess) > 0
        max_guess = guess;
        if max_guess==1
            res = 0;
            return;
        end
    elseif results(guess) < 0
        min_guess = guess;
    end
    
    if min_guess+1 >= max_guess
        res = mu_test(min_guess);
        return
    end
    count = count+1;
end
error('Not sure how we got here');
end

function res = test_mu(mu, dt, M2, M3)
rng(0);
global Z
M2plus = M2;
M3plus = M3;
M2minus = -M2;
M3minus = -M3;
N_max = 500;

Je = inv(Z(1:3,1:3));
w_max = pi/180; % 1 deg/s absolute limit on the largest axis
E_max = [1,0,0]*Je*[1;0;0]*w_max^2;
wheel_rate_limit = 2*pi*6000/60; % the max wheel rate in radians per second, note this leads to poor numerical conditioning
w_max = pi/180; % 1 deg/s absolute limit on the largest axis
year = 366*24*3600; % this is just used for scaling the time variable

p_con = @(x) rh(x(1)*year, x(2:7)); % in the safe set of quaternions
w_con = @(x) (x(5:7)'*Je*x(5:7) - E_max)*1e3; % in the safe set of angular velocities
r_con = @(x) norm(x(2:4)) - 1; 

%% Find delta_1 when symmetric
delta1 = 1/2*(mu+M2plus-M2minus)*(dt/2)^2 - 1/6*M3minus*(dt/2)^3;
    % this formula is specific to the symmetric case (M2plus = M2minus, etc)

%% Find Delta3
N = 100;
f_left = @(gamma, tau1) -tau1^3*((M3minus*(M2plus - M2minus + mu))/(3*mu) + (M3plus*(M2minus - M2plus))/(6*mu)) - tau1^2*(((M2minus - M2plus)*(M2plus - M2minus + mu))/(2*mu) + (M3minus*gamma)/(2*mu)) - (gamma*tau1*(M2minus - M2plus))/mu - (M3minus*M3plus*tau1^4)/(8*mu);
f_right = @(gamma, tau1, tau2) ((M3minus*((M2minus - M2plus)/mu - 2))/6 - (M3minus*(M2plus - M2minus + mu))/(3*mu))*(dt - tau1)^3 - tau2^3*((M3minus*((M2minus - M2plus)/mu - 2))/6 - (M3minus*(M2plus - M2minus + mu))/(3*mu)) + tau2^2*((((M2minus - M2plus)/mu - 2)*(M2plus - M2minus + mu))/2 - (M3minus*gamma)/(2*mu)) - ((((M2minus - M2plus)/mu - 2)*(M2plus - M2minus + mu))/2 - (M3minus*gamma)/(2*mu))*(dt - tau1)^2 + tau2^3*((M3minus*(M2minus - M2plus + mu))/(3*mu) - (M3plus*(M2minus - M2plus))/(6*mu)) + tau2^2*(((M2minus - M2plus)*(M2minus - M2plus + mu))/(2*mu) - (M3minus*gamma)/(2*mu)) - gamma*tau2*((M2minus - M2plus)/mu - 2) + gamma*(dt - tau1)*((M2minus - M2plus)/mu - 2) - (M3minus^2*tau2^4)/(8*mu) + (M3minus^2*(dt - tau1)^4)/(8*mu) - (gamma*tau2*(M2minus - M2plus))/mu - (M3minus*M3plus*tau2^4)/(8*mu); 

scale = 1e4;
func = @(x) -min(f_left(x(1)/scale, x(2)), f_right(x(1)/scale, x(2), x(3)))*1e6;
lower = [0; 0; 0];
M1max = (-mu-M2plus+M2minus)*(-dt) + M3minus*(-dt^2)/2;
upper = [M1max; dt; dt];
rule1 = @(x) x(3) - ( (-mu+M2plus-M2minus)^2 - 2*M3plus*x(1) >= 0 )*(-(-mu+M2plus-M2minus)-sqrt( (-mu+M2plus-M2minus)^2 - 2*M3plus*x(1) ))/M3plus ...
    - ( (-mu+M2plus-M2minus)^2 - 2*M3plus*x(1) < 0 )*(dt - x(2));
rule2 = @(x) -x(3) + (-(-mu+M2minus-M2plus)-sqrt( (-mu+M2minus-M2plus)^2 - 2*M3minus*x(1) ))/M3minus;
nonlcon = @(x) assemble_constraint(x, @(y) [rule1(y./[scale;1;1]); rule2(y./[scale;1;1])]);
min_val = 0;
min_res = [];
upper(1) = upper(1)*scale;
A = [0, 1, 1];
for i=1:min(N, N_max)
guess = rand(3,1).*upper;
[res, fval] = fmincon(func, guess, A, dt, [], [], lower, upper, nonlcon, optimset('Display', 'off'));
if A*res-dt<1e-5 && max(nonlcon(res)) < 1e-5
    min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
end
min_val = min_val/1e6;

Delta3 = -min_val;

%% Fast Method for checking Delta_2, delta_2
% delta_2 = 0.97e-5; % This is an asymmetric case that works, but provides minimal benefit
% Delta_2 = 1.3e-5;
my_ceil = @(x, n) ceil(x*(10^n))/((10)^n);
delta_2 = my_ceil(delta1, 8);
Delta_2 = max(my_ceil(delta1, 8), my_ceil(Delta3, 8));

%% Look at the future
h_future = @(x) rh(x(1)*year, x(2:7)) + rhdot(x(1)*year, x(2:7))*dt ...
    + 1/2*rphi(x(1)*year, [x(2:7); x(8:11)*1e3], x(12:15))*dt^2 ...
    + 1/2*M2plus*dt^2 + 1/6*M3plus*dt^3;
hdot_future = @(x) rhdot(x(1)*year, x(2:7)) ...
    + rphi(x(1)*year, [x(2:7); x(8:11)*1e3], x(12:15))*dt ...
    + M2plus*dt + 1/2*M3plus*dt^2;
H_future = @(x) h_future(x) + 1/(2*mu)*absSq(hdot_future(x));

wheel_limit = 0.001; % https://www.cubespace.co.za/products/adcs-components/cubewheel/#cubewheel-specifications
wheel_limit = 0.7*wheel_limit;
u_guess = @(x) sign(get_s(x(1)*year)'*skew(x(2:4))*Z(1:3,4:7))'*wheel_limit;
h_future_star = @(x) h_future([x(1:11); u_guess(x)]);
H_future_star = @(x) H_future([x(1:11); u_guess(x)]);

func1 = @(x) -h_future_star(x)*1e4;
func2 = @(x) -H_future_star(x)*1e4;
lower = -[0; ones(3,1); 3*w_max*ones(3,1); wheel_rate_limit*ones(4,1)/1e3];
upper =  [1; ones(3,1); 3*w_max*ones(3,1); wheel_rate_limit*ones(4,1)/1e3];
H_con = @(x) rh(x(1)*year,x(2:7)) + 1/(2*mu)*absSq(rhdot(x(1)*year,x(2:7))) + Delta_2;
nonlcon = @(x) assemble_constraint(x, @(y) [(p_con(y)+delta_2)*1e3; H_con(y)*1e4; w_con(y)*1e3], @(y) r_con(y));
N = 150;
min_val1 = Inf;
min_res1 = [];
min_val2 = Inf;
min_res2 = [];
N_pass1 = 0;
N_pass2 = 0;
options = optimset('Display', 'off', 'OutputFcn', @fmincon_stop);
for i=1:min(N, N_max)
guess = randn(11,1).*upper;
tic;
[res, fval] = fmincon(func1, guess, [], [], [], [], lower, upper, nonlcon, options);
if sum(nonlcon(res) > 0) == 0, min_val1 = min(min_val1, fval); if min_val1 == fval, min_res1 = res; end; N_pass1 = N_pass1+1; end
tic;
[res, fval] = fmincon(func2, guess, [], [], [], [], lower, upper, nonlcon, options); waitbar(i/N)
if sum(nonlcon(res) > 0) == 0, min_val2 = min(min_val2, fval); if min_val2 == fval, min_res2 = res; end; N_pass2 = N_pass2+1; end
end
min_val1 = -min_val1/1e4 + delta_2;
min_val2 = -min_val2/1e4 + Delta_2;

% If both of the above are negative, then we found a valid mu
% If both of the above are much less than zero, then we can increase mu
% If min_val2 is much less than zero while min_val1 is not, then we can decrease delta_2
% (and potentially increase Delta_2 at the same time) which may potentially allow us to 
% increase mu. Higher mu means bigger safe set.
% An "optimally" tuned set of constants delta_2, Delta_2, mu will be such that delta_2 is
% decreased from delta1 and Delta_2 is increased from delta1 until min_val1 is
% approximately equal to min_val2. Then we increase mu until min_val1 and min_val2 are
% both zero.
% If min_val2 is less than min_val1 when Delta_2 = delta1, then you may want to rerun this
% section, because that does not often happen.

res = max(min_val1, min_val2);

% 0.00167 is the largest mu to three significant digits for which we could find any
% combination of delta_2 and Delta_2 that was valid.

end

function stop = fmincon_stop(x,optimValues,state)
stop = toc > 10;
if stop
    disp('fmincon force stopped');
end
end
%%
function u = u_verify(t, x, constants)
x(8:11) = x(8:11)*1e3; % adjust wheel velocities to match what this function uses

% Relative Degree 2 Constraint
[A1, b1] = ConstraintQ(t, x, constants, constants.p1);
[A2, b2] = ConstraintQ(t, x, constants, constants.p2);

% Relative Degree 1 Constraint
[H0, A0, b0] = ConstraintE(t, x, constants);
% The above constraints will output warnings that can be ignored in this case.

lower = -constants.wheel_limit*[1;1;1;1];
upper =  constants.wheel_limit*[1;1;1;1];

u = linprog_exact(A0, [A1; A2], [b1; b2], lower, upper);
end

function out = h_feas(t,x,p)
global O_Earth_from_Sol omega_Earth ctheta
r = p;
s = QtoR(x(1:4))'*O_Earth_from_Sol*[cos(omega_Earth*t); sin(omega_Earth*t); 0];

out = dot(s,r) - ctheta;
end

function out = hdot_feas(t,x,p)
global O_Earth_from_Sol omega_Earth

s = QtoR(x(1:4))'*O_Earth_from_Sol*[cos(omega_Earth*t); sin(omega_Earth*t); 0];
sdot = QtoR(x(1:4))'*O_Earth_from_Sol*[-omega_Earth*sin(omega_Earth*t); omega_Earth*cos(omega_Earth*t); 0];
r = p;
w = x(5:7);

out = sdot'*r + s'*skew(w)*r;
end

function out = phi_feas(t,x,u,p)
global O_Earth_from_Sol omega_Earth Z Jtot wheel_axis Jw
s = QtoR(x(1:4))'*O_Earth_from_Sol*[cos(omega_Earth*t); sin(omega_Earth*t); 0];
sdot = QtoR(x(1:4))'*O_Earth_from_Sol*[-omega_Earth*sin(omega_Earth*t); omega_Earth*cos(omega_Earth*t); 0];
sddot = QtoR(x(1:4))'*O_Earth_from_Sol*[-omega_Earth^2*cos(omega_Earth*t); -omega_Earth^2*sin(omega_Earth*t); 0];
r = p;
w = x(5:7);
W = x(8:11);
wdot = Z(1:3,:)*[-cross(w,Jtot*w + wheel_axis*Jw*W); u];

out = sddot'*r + 2*sdot'*skew(w)*r + s'*skew(w)^2*r - s'*skew(r)*wdot;
end