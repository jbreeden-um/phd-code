%% Mass Parameter Setup
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
x = fsolve(f, [0; 0.5; 0.5]);
tetrahedron_points(4,:) = x';
global wheel_axis
wheel_axis(:,1) = cross(tetrahedron_points(1,:) - tetrahedron_points(2,:), tetrahedron_points(1,:) - tetrahedron_points(3,:))';
wheel_axis(:,2) = cross(tetrahedron_points(1,:) - tetrahedron_points(4,:), tetrahedron_points(1,:) - tetrahedron_points(2,:))';
wheel_axis(:,3) = cross(tetrahedron_points(1,:) - tetrahedron_points(3,:), tetrahedron_points(1,:) - tetrahedron_points(4,:))';
wheel_axis(:,4) = cross(tetrahedron_points(3,:) - tetrahedron_points(4,:), tetrahedron_points(2,:) - tetrahedron_points(3,:))';
    % the above has been ordered so that all the wheels positive spin axes face away from 
    % the center of the spacecraft
wheel_axis = wheel_axis./vecnorm(wheel_axis)
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
torque = force*lever;

global O_Earth_from_Sol omega_Earth
O_Earth_from_Sol = [1, 0, 0; 0, cosd(-23.4392911), sind(-23.4392911); 0, -sind(-23.4392911), cosd(-23.4392911)]; % Ref: https://github.com/ericstoneking/42/blob/abd2b0bfbee8c10bf86e07d726340ef0bf0f1308/Source/42init.c
    % Orientation matrix converts vectors from ecliptic frame to Earth polar frame
omega_Earth = 2*pi/(365.25*24*3600);

dt = 0.2;
year = 366*24*3600; % this is just used for scaling the time variable
N_max = 500;

%% 
N_max = 0; % this line only runs when you run the whole program through rather than one section at a time.

%% Find M_2
tic
g1 = [zeros(3,4); zeros(3,4); eye(4)];
g2 = [zeros(3,3); Z(1:3,1:3); zeros(4,3)];

% x(1) = t, x(2:4) = p, x(5:7) = omega, x(8:10) = d
func = @(x) rgrad_hdot(x(1)*year, x(2:7))*g2*x(8:10);
p_con = @(x) rh(x(1)*year, x(2:7)); % in the safe set of quaternions
w_con = @(x) (x(5:7)'*Je*x(5:7) - E_max)*1e3; % in the safe set of angular velocities
d_con = @(x) norm(x(8:10)) - torque; % disturbance within limit
r_con = @(x) norm(x(2:4)) - 1; 
    % This should be an equality constraint, but the code is more efficient as is.
    % Note that later on we will enforce this equality, but here the worst case occurs
    % with equality anyways.
lower = -[0; ones(3,1); 3*w_max*ones(3,1); torque*ones(3,1)];
upper =  [1; ones(3,1); 3*w_max*ones(3,1); torque*ones(3,1)];
nonlcon = @(x) assemble_constraint(x, @(y) [p_con(y); w_con(y); d_con(y); r_con(y)]);
N = 10;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(10,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val
if ~isempty(min_res), nonlcon(min_res), end

M2minus = -Z(3,3)*torque;
toc

%%
func = @(x) -rgrad_hdot(x(1)*year, x(2:7))*g2*x(8:10);
N = 10;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(10,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val
if ~isempty(min_res), nonlcon(min_res), end

M2plus = Z(3,3)*torque;
M2 = max(M2plus, abs(M2minus));

%% Find M_3
tic
% omega_Earth = 4e-6; 
    % This is to test what happens if the time-varying component is sped up, such as would
    % occur for a lunar constraint, or if we include orbital dynamics about the Earth
    % instead of just about the Sun.
    % This rate is still so small as to make negligible difference in the value of M3.

% x(1) = t, x(2:4) = p, x(5:7) = omega, x(8:10) = d, x(11:14) = W, x(15:18) = wheel_commands
func = @(x) rphidot(x(1)*year, [x(2:7); x(11:14)*1e3], x(15:18), x(8:10));
lower = -[0; ones(3,1); w_max*ones(3,1)*3; torque*ones(3,1); wheel_rate_limit*ones(4,1)/1e3; wheel_limit*ones(4,1)];
upper =  [1; ones(3,1); w_max*ones(3,1)*3; torque*ones(3,1); wheel_rate_limit*ones(4,1)/1e3; wheel_limit*ones(4,1)];
nonlcon = @(x) assemble_constraint(x, @(y) [p_con(y); w_con(y); d_con(y); r_con(y)]);
N = 100;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(18,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val
if ~isempty(min_res), nonlcon(min_res), end

M3minus = -0.0062;
toc

%%

func = @(x) -rphidot(x(1)*year, [x(2:7); x(11:14)*1e3], x(15:18), x(8:10));
N = 100;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(18,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val
if ~isempty(min_res), nonlcon(min_res), end

M3plus = 0.0062;
M3 = max(M3plus, abs(M3minus));

%% Comparison to naive overshoot approach
% Second derivative of H with no noise
muC = 0.00167;
func = @(x) ( rphi(x(1)*year, [x(2:7); x(11:14)*1e3], x(15:18)) ...
    + rphi(x(1)*year, [x(2:7); x(11:14)*1e3], x(15:18))^2*sign(rhdot(x(1)*year, x(2:7)))/muC ...
    + abs(rhdot(x(1)*year, x(2:7)))*rphidot(x(1)*year, [x(2:7); x(11:14)*1e3], x(15:18), [0;0;0])/muC );
lower = -[0; ones(3,1); w_max*ones(3,1)*3; torque*ones(3,1); wheel_rate_limit*ones(4,1)/1e3; wheel_limit*ones(4,1)];
upper =  [1; ones(3,1); w_max*ones(3,1)*3; torque*ones(3,1); wheel_rate_limit*ones(4,1)/1e3; wheel_limit*ones(4,1)];
nonlcon = @(x) assemble_constraint(x, @(y) [p_con(y); w_con(y); d_con(y); r_con(y)]);
N = 100;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(18,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
if sum(nonlcon(res) > 0)==0, min_val = min(min_val, fval); if min_val == fval, min_res = res; end; end
end
min_val
if ~isempty(min_res), nonlcon(min_res), end

Hddot_min = -0.5504;
Delta = Hddot_min*dt^2/8;

%%
% Second derivative of h with no noise
u_guess = @(x) sign(s_func(x(1)*year)'*skew(x(2:4))*Z(1:3,4:7))'*wheel_limit; % control law that most decrease h
lower = -[0; ones(3,1); w_max*ones(3,1)*3; wheel_rate_limit*ones(4,1)/1e3];
upper =  [1; ones(3,1); w_max*ones(3,1)*3; wheel_rate_limit*ones(4,1)/1e3];
func = @(x) rphi(x(1)*year, [x(2:7); x(8:11)*1e3], u_guess(x))*1e6;
N = 100;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(11,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val = min_val / 1e6
if ~isempty(min_res), nonlcon(min_res), end

hddot_min = -0.0262;
Delta = hddot_min*dt^2/8;

%%
% Third Derivative of h with no noise
func = @(x) rphidot(x(1)*year, [x(2:7); x(11:14)*1e3], x(15:18), [0;0;0]);
N = 100;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(18,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val
if ~isempty(min_res), nonlcon(min_res), end

hdddot_min = -0.0062; % the disturbance has a negligible effect on this constant, unlike M2plus/M2minus
hddot_eff = muC - hdddot_min*dt;
Delta = hddot_eff*dt^2/8;

%% Comparison to continuous time
muC = 0.0025;
u_guess = @(x) sign(get_s(x(1)*year)'*skew(x(2:4))*Z(1:3,4:7))'*wheel_limit;
H_conC = @(x) rh(x(1)*year, x(2:7)) + 1/(2*muC)*absSq(rhdot(x(1)*year, x(2:7)));
Hdot = @(x) rhdot(x(1)*year, x(2:7)) + 1/muC*abs(rhdot(x(1)*year, x(2:7))) ...
    *rhddot(x(1)*year, [x(2:7); x(11:14)*1e3], u_guess(x), x(8:10));

func = @(x) -Hdot(x)*1e6;
lower = -[0; ones(3,1); w_max*ones(3,1)*3; ones(3,1)*torque; wheel_rate_limit*ones(4,1)/1e3];
upper =  [1; ones(3,1); w_max*ones(3,1)*3; ones(3,1)*torque; wheel_rate_limit*ones(4,1)/1e3];
nonlcon = @(x) assemble_constraint(x, @(y) [p_con(y); w_con(y)*1e3; d_con(y)*1e3], @(y) [r_con(y)*1e4; H_conC(y)*1e2]);
N = 200;
min_val = inf;
min_res = [];
for i=1:min(N, N_max)
if i==0, guess = saved; else, guess = randn(14,1).*upper; end
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
if sum(nonlcon(res) > 0) == 0
    min_val = min(min_val, fval); if min_val == fval, min_res = res; end
    if min_val < 0, break; end
end
end
min_val = min_val / 1e6
if ~isempty(min_res), nonlcon(min_res), end

% The selection of muC above is only accurate to 2 significant digits.
% That is, 0.0025 worked, and 0.0026 did not work.
% To test this yourself, put in a guess for muC, and run this section. 
% If min_val is positive, then muC is accepted. If min_val is negative, then muC is too
% large. You may want to run this code multiple times to be confident muC is viable.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following segment of this file uses the above constants to determine a valid set of
% mu, delta_2, Delta_2 constants.

%% Start finding Delta, delta, and mu
mu = 0.00167; % our guess for mu

%% Find delta_1 when symmetric
delta1 = 1/2*(mu+M2plus-M2minus)*(dt/2)^2 - 1/6*M3minus*(dt/2)^3
    % this formula is specific to the symmetric case (M2plus = M2minus, etc)

%% Find delta_1 when assymetric
kappa2_guess = @(sigma) 1/2*(mu+M2plus-M2minus)*(dt-sigma)^2 - 1/6*M3minus*(dt-sigma)^3;
kappa1_guess = @(sigma) 1/2*(mu+M2plus-M2minus)*(sigma)^2 + 1/6*M3plus*(sigma)^3;
common_val = fzero(@(s) kappa2_guess(s)-kappa1_guess(s), 0.75*dt);
delta1_check = kappa2_guess(common_val)

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
[res, fval] = fmincon(func, guess, A, dt, [], [], lower, upper, nonlcon);
if A*res-dt<1e-5 && max(nonlcon(res)) < 1e-5
    min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
end
min_val = min_val/1e6

Delta3 = 1.09e-5;

%% Fast Method for checking Delta_2, delta_2
% delta_2 = 0.97e-5; % This is an asymmetric case that works, but provides minimal benefit
% Delta_2 = 1.3e-5;
my_ceil = @(x, n) ceil(x*(10^n))/((10)^n);
delta_2 = my_ceil(delta1, 8);
Delta_2 = my_ceil(delta1, 8);
kappa2 = @(sigma) -delta_2 + 1/2*(mu+M2plus-M2minus)*(dt-sigma)^2 - 1/6*M3minus*(dt-sigma)^3;
kappa1 = @(sigma) -Delta_2 + 1/2*(mu+M2plus-M2minus)*(sigma)^2 + 1/6*M3plus*(sigma)^3;

common_val = fzero(@(s) kappa2(s)-kappa1(s), 0.75*dt);
check = kappa2(common_val)
% If check <= 0, then this combination of delta and Delta_b is accepted
% The above method usually works, but when delta_2 and Delta_2 are very different, it can
% become inconsistent.

%% Complete Method for checking Delta_2, delta_2
func = @(x) max(kappa1(x), kappa2(x))*1e5;
min_val = inf;
min_res = [];
lower = 0;
upper = dt;
N = 10;
for i=1:min(N, N_max)
guess = rand*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
check = min_val / 1e5

% If check <= 0, then delta_2 and Delta_2 are accepted.
% Conservatism is reduced when check is closest to zero
% For this problem, 1e-9 to 1-8 are considered "close" to zero.

if Delta_2 < Delta3
    warning('Results may be inaccurate because Delta2 tolerance is not large enough');
end

%% Look at the future
h_future = @(x) rh(x(1)*year, x(2:7)) + rhdot(x(1)*year, x(2:7))*dt ...
    + 1/2*rphi(x(1)*year, [x(2:7); x(8:11)*1e3], x(12:15))*dt^2 ...
    + 1/2*M2plus*dt^2 + 1/6*M3plus*dt^3;
hdot_future = @(x) rhdot(x(1)*year, x(2:7)) ...
    + rphi(x(1)*year, [x(2:7); x(8:11)*1e3], x(12:15))*dt ...
    + M2plus*dt + 1/2*M3plus*dt^2;
H_future = @(x) h_future(x) + 1/(2*mu)*absSq(hdot_future(x));

h_future_star = @(x) h_future([x(1:11); u_guess(x)]);
H_future_star = @(x) H_future([x(1:11); u_guess(x)]);

func1 = @(x) -h_future_star(x)*1e4;
func2 = @(x) -H_future_star(x)*1e4;
lower = -[0; ones(3,1); 3*w_max*ones(3,1); wheel_rate_limit*ones(4,1)/1e3];
upper =  [1; ones(3,1); 3*w_max*ones(3,1); wheel_rate_limit*ones(4,1)/1e3];
H_con = @(x) rh(x(1)*year,x(2:7)) + 1/(2*mu)*absSq(rhdot(x(1)*year,x(2:7))) + Delta_2;
nonlcon = @(x) assemble_constraint(x, @(y) [(p_con(y)+delta_2)*1e3; H_con(y)*1e4; w_con(y)*1e3], @(y) r_con(y));
N = 300;
min_val1 = Inf;
min_res1 = [];
min_val2 = Inf;
min_res2 = [];
N_pass1 = 0;
N_pass2 = 0;
for i=1:min(N, N_max)
guess = randn(11,1).*upper;
[res, fval] = fmincon(func1, guess, [], [], [], [], lower, upper, nonlcon);
if sum(nonlcon(res) > 0) == 0, min_val1 = min(min_val1, fval); if min_val1 == fval, min_res1 = res; end; N_pass1 = N_pass1+1; end
[res, fval] = fmincon(func2, guess, [], [], [], [], lower, upper, nonlcon);
if sum(nonlcon(res) > 0) == 0, min_val2 = min(min_val2, fval); if min_val2 == fval, min_res2 = res; end; N_pass2 = N_pass2+1; end
end
min_val1 = -min_val1/1e4 + delta_2
min_val2 = -min_val2/1e4 + Delta_2

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

% 0.00167 is the largest mu to three significant digits for which we could find any
% combination of delta_2 and Delta_2 that was valid.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The remainder of this file compares the effective margins for both types of h_omega functions

%% Two Norm Rate Constraint WITHOUT Compensating for Controllable Terms of PhiDot
Z11 = Z(1:3,1:3);
Z12 = Z(1:3, 4:7);
Z21 = Z(4:7,1:3);
Z22 = Z(4:7,4:7);
wdot = @(w, W, u) Z11*(-skew(w)*(Jtot*w+wheel_axis*Jw*W)) + Z12*u;
wdot_v = @(w, W, u, v) Z11*(v-skew(w)*(Jtot*w+wheel_axis*Jw*W)) + Z12*u;
Wdot_v = @(w, W, u, v) Z21*(v-skew(w)*(Jtot*w+wheel_axis*Jw*W)) + Z22*u;
phidot_twonorm = @(w, W, u, v) 2*wdot(w,W,u)'*wdot_v(w,W,u,v) ...
    + 2*w'*(Z11*(-skew(wdot_v(w,W,u,v))*(Jtot*w+wheel_axis*Jw*W) ...
        - skew(w)*(Jtot*wdot_v(w,W,u,v)+wheel_axis*Jw*Wdot_v(w,W,u,v))));

func = @(x) -phidot_twonorm(x(1:3), x(11:14)*1e3, x(7:10), x(4:6))*1e3;
% Note that the following uses applies w_max to all three axes, and thus induces a smaller
% safe set then the energy constraint method. We could also have applied 1.65*w_max to all
% axes to create a larger safe set.
lower = -[w_max*ones(3,1); torque*ones(3,1); wheel_limit*ones(4,1); wheel_rate_limit/1e3*ones(4,1)];
upper =  [w_max*ones(3,1); torque*ones(3,1); wheel_limit*ones(4,1); wheel_rate_limit/1e3*ones(4,1)];
w_con = @(x) norm(x(1:3)) - w_max;
d_con = @(x) norm(x(4:6)) - torque;
nonlcon = @(x) assemble_constraint(x, @(y) [w_con(y)*1e3; d_con(y)*1e3]); % this line helps with tolerances but does not seem to affect the final result
N = 100;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(14,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val = min_val / 1e3
if ~isempty(min_res), nonlcon(min_res), end

percent_margin0 = min_val*dt^2/2 / w_max^2
% 10% margin is pretty bad
% Note the proportion of conservatism is smaller when w_max is larger.

%% Two Norm Rate Constraint WITH Compensating for Controllable Term of PhiDot
error_twonorm = @(x) phidot_twonorm(x(1:3), x(11:14)*1e3, x(7:10), x(4:6)) - 2*x(7:10)'*Z12'*Z12*x(7:10);

func = @(x) -error_twonorm(x)*1e3;
N = 100;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(14,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val = min_val / 1e3
if ~isempty(min_res), nonlcon(min_res), end

percent_margin1 = min_val*dt^2/2 / w_max^2
% 3.2% margin is acceptable, but not great.

%% Energy Constraint WITH Compensating for Controllable Term of PhiDot
phidot_energy = @(w, W, u, v) 2*(Z11*(v-skew(w)*(Jtot*w+wheel_axis*Jw*W))+Z12*u)'*inv(Z11)*Z12*u;
error_energy = @(x) phidot_energy(x(1:3), x(11:14)*1e3, x(7:10), x(4:6)) - 2*x(7:10)'*Z12'*inv(Z11)*Z12*x(7:10);

func = @(x) -error_energy(x)*1e4;
lower = -[3*w_max*ones(3,1); torque*ones(3,1); wheel_limit*ones(4,1); wheel_rate_limit/1e3*ones(4,1)];
upper =  [3*w_max*ones(3,1); torque*ones(3,1); wheel_limit*ones(4,1); wheel_rate_limit/1e3*ones(4,1)];
w_con = @(x) x(1:3)'*Je*x(1:3) - E_max; 
nonlcon = @(x) assemble_constraint(x, @(y) [w_con(y)*1e3; d_con(y)*1e3]);
N = 100;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(14,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val = min_val/1e4
if ~isempty(min_res), nonlcon(min_res), end

percent_margin2 = min_val*dt^2/2 / E_max
% 0.77% margin is very good.

M2omega = 1.95e-5;

%% Energy Constraint WITHOUT Compensating for Controllable terms of PhiDot
func = @(x) -phidot_energy(x(1:3), x(11:14)*1e3, x(7:10), x(4:6))*1e5; 
N = 100;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(14,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val = min_val/1e5
if ~isempty(min_res), nonlcon(min_res), end

percent_margin3 = min_val*dt^2/2 / E_max
% 3.3% margin is okay.
% If we wish to keep the problem simply a QP and not a QCQP, then this is much better than
% the 10% margin introduced by the two norm rate constraint

M2omega_lin = 8.30e-5;

%% Applying LCSS Method 2 To Two Norm Rate Constraint
% This and the following section is purely for comparison to a method not actually
% presented in this paper, but which can be easily derived from the work in our prior
% L-CSS paper (see \phi_2^g and \nu_2^g).
% Note that this method introduces a more complex optimization problem to compute the 
% necessary constants, which has been troublesome on certain systems.
psi_omega = @(w, W, u, v) 2*w'*(Z11*(v-skew(w)*(Jtot*w+wheel_axis*Jw*W)) + Z12*u);

% x(1:3) = w1, x(4:6) = w2, x(7:9) = d1, x(10:12) = d2, x(13:16) = wheels
difference = @(x) psi_omega(x(1:3), x(7:10)*1e3, x(15:18), x(19:21)) ...
    - psi_omega(x(4:6), x(11:14)*1e3, x(15:18), x(22:24));
% Note that this function omits the class-k function used in the L-CSS paper. One can
% prove that it is unnecessary using the same logic as used in this paper.

func = @(x) -norm(wdot_v(x(1:3), x(11:14)*1e3, x(7:10), x(4:6)))*1e3;
lower = -[w_max*ones(3,1); torque*ones(3,1); wheel_limit*ones(4,1); wheel_rate_limit/1e3*ones(4,1)];
upper =  [w_max*ones(3,1); torque*ones(3,1); wheel_limit*ones(4,1); wheel_rate_limit/1e3*ones(4,1)];
w_con = @(x) norm(x(1:3)) - w_max;
d_con = @(x) norm(x(4:6)) - torque;
nonlcon = @(x) assemble_constraint(x, @(y) [w_con(y)*1e3; d_con(y)*1e3]);
[~, min_val] = fmincon(func, zeros(14,1), [], [], [], [], lower, upper, nonlcon);
max_w_change = (-min_val/1e3)*dt;
% This number is quite large, which partially explains why this method yields such poor
% results. Maybe we need to propagate reachable sets more directly.

func = @(x) -norm(Wdot_v(x(1:3), x(11:14)*1e3, x(7:10), x(4:6)))*1e3;
[~, min_val] = fmincon(func, [0;1;0;0;0;1;-1;1;1;1;1;1;-1;1].*upper, [], [], [], [], lower, upper, nonlcon);
max_W_change = (-min_val/1e3)*dt;

func = @(x) difference(x)*1e4;
lower = -[w_max*ones(6,1); wheel_rate_limit/1e3*ones(8,1); wheel_limit*ones(4,1); torque*ones(6,1)];
upper =  [w_max*ones(6,1); wheel_rate_limit/1e3*ones(8,1); wheel_limit*ones(4,1); torque*ones(6,1)];
w_con1 = @(x) norm(x(1:3)) - w_max;
w_con2 = @(x) norm(x(4:6)) - w_max;
d_con1 = @(x) norm(x(19:21)) - torque;
d_con2 = @(x) norm(x(22:24)) - torque;
reach_con1 = @(x) norm(x(1:3)-x(4:6)) - max_w_change;
reach_con2 = @(x) norm(x(7:10)-x(11:14))*1e3 - max_W_change;
nonlcon = @(x) assemble_constraint(x, @(y) [w_con1(y); w_con2(y); d_con1(y); d_con2(y); reach_con1(y); reach_con2(y)/1e6]*1e3);
N = 100;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(24,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val = min_val / 1e4
if ~isempty(min_res), nonlcon(min_res), end

percent_margin = min_val*dt / w_max^2
% 20.5% margin is pretty bad.

%% Applying LCSS Method 2 To Energy Constraint
psi_omega = @(w, u, v) 2*w'*v + 2*w'*inv(Z11)*Z12*u; % this derivative is not a function of W, yay!

% x(1:3) = w1, x(4:6) = w2, x(7:9) = d1, x(10:12) = d2, x(13:16) = wheels
difference = @(x) psi_omega(x(1:3), x(13:16), x(7:9)) - psi_omega(x(4:6), x(13:16), x(10:12));

func = @(x) difference(x)*1e4;
lower = -[3*w_max*ones(6,1); torque*ones(6,1); wheel_limit*ones(4,1)];
upper =  [3*w_max*ones(6,1); torque*ones(6,1); wheel_limit*ones(4,1)];
w_con = @(x) [x(1:3)'*Je*x(1:3); x(4:6)'*Je*x(4:6)] - E_max; 
d_con = @(x) [norm(x(7:9)); norm(x(10:12))] - torque;
reach_con = @(x) norm(x(1:3)-x(4:6)) - max_w_change;
nonlcon = @(x) assemble_constraint(x, @(y) [w_con(y); d_con(y); reach_con(y)]*1e3);
N = 100;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(16,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val = min_val / 1e4
if ~isempty(min_res), nonlcon(min_res), end

percent_margin = min_val*dt / E_max
% 7.8% margin is better, but still poor compared to the higher order methods, even without
% using the additional QCQP term.

%% Find M_1omega
func = @(x) -2*x(1:3)'*x(4:6)*1e8;
w_con = @(x) x(1:3)'*Je*x(1:3) - E_max; 
d_con = @(x) norm(x(4:6)) - torque;
lower = -[3*w_max*ones(3,1); torque*ones(3,1)];
upper =  [3*w_max*ones(3,1); torque*ones(3,1)];
nonlcon = @(x) assemble_constraint(x, @(y) [w_con(y); d_con(y)]*1e6);
N = 10;
min_val = 0;
min_res = [];
for i=1:min(N, N_max)
guess = randn(6,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val = min_val / 1e8
if ~isempty(min_res), nonlcon(min_res), end

M1omega = 2*w_max*torque*sqrt(Je(1,1)/Je(3,3));

%% Check for feasibility of Energy Constraint
h_omega = @(x) x(1:3)'*Je*x(1:3) - E_max;
psi_omega = @(w, u) 2*w'*inv(Z11)*Z12*u;
psidot_omega = @(w, W, u, v) 2*(Z11*(v-skew(w)*(Jtot*w+wheel_axis*Jw*W))+Z12*u)'*inv(Z11)*Z12*u;
psidot_omega_known = @(u) 2*u'*Z12'*inv(Z11)*Z12*u;
u_omega_opt = @(w) -sign(w'*inv(Z11)*Z12)'*wheel_limit;
h_future = @(x) h_omega(x) + (psi_omega(x(1:3),u_omega_opt(x(1:3))) + M1omega)*dt ...
    + (psidot_omega_known(u_omega_opt(x(1:3))) + M2omega)*dt^2/2;
    % This only works because psidot_omega_known is positive definite!
func = @(x) -h_future(x)*1e4;
lower = -[3*w_max*ones(3,1)];
upper =  [3*w_max*ones(3,1)];
nonlcon = @(x) assemble_constraint(x, @(y) [w_con(y)]*1e6);
N = 100;
min_val = 1;
min_res = [];
for i=1:min(N, N_max)
guess = randn(3,1).*upper;
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
min_val = min_val/1e4
if ~isempty(min_res), nonlcon(min_res), end

% min_val was 6.47e-6.
% If min_val is positive, then feasibility is guaranteed. If min_val is not positive, then
% something is probably broekn in our problem formulation (e.g. dt is too large).

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The next section takes a while to run and does not produce anything required parameters,
% so I'm cutting off the file here. Feel free to run the next section manually.
%% Check for joint feasibility of all the constraints
% Ignore all the Excessive H, excessive kinetic energy warnings that appear. This is just
% because each optimization may start outside the safe set.
constants = load('parameters.mat');
state = @(x) [x(1:4); x(5:7); x(8:11)*1e3];

p_con1 = @(x) h_feas(x(1)*year, state(x(2:12)), constants.p1) + delta_2;
p_con2 = @(x) h_feas(x(1)*year, state(x(2:12)), constants.p2) + delta_2;
w_con = @(x) (x(6:8)'*Je*x(6:8) - E_max);
r_con = @(x) norm(x(2:5)) - 1;
lower = -[0; ones(4,1); 3*w_max*ones(3,1); wheel_rate_limit/1e3*ones(4,1)];
upper =  [1; ones(4,1); 3*w_max*ones(3,1); wheel_rate_limit/1e3*ones(4,1)];
H_con1 = @(x) h_feas(x(1)*year,state(x(2:12)),constants.p1) + 1/(2*mu)*absSq(hdot_feas(x(1)*year,state(x(2:12)),constants.p1)) + Delta_2;
H_con2 = @(x) h_feas(x(1)*year,state(x(2:12)),constants.p2) + 1/(2*mu)*absSq(hdot_feas(x(1)*year,state(x(2:12)),constants.p2)) + Delta_2;
nonlcon = @(x) assemble_constraint(x, ...
    @(y) [[p_con1(y); p_con2(y)]*1e1; [H_con1(y); H_con2(y)]*1e1; w_con(y)*1e6], @(y) r_con(y));

h_omega_future = @(x, u) h_omega(x(6:8)) + (psi_omega(x(6:8),u) + M1omega)*dt ...
    + (psidot_omega_known(u) + M2omega)*dt^2/2;
h_future1 = @(x, u) h_feas(x(1)*year, state(x(2:12)), constants.p1) ...
    + hdot_feas(x(1)*year, state(x(2:12)), constants.p1)*dt ...
    + 1/2*phi_feas(x(1)*year, state(x(2:12)), u, constants.p1)*dt^2 ...
    + 1/2*M2plus*dt^2 + 1/6*M3plus*dt^3;
hdot_future1 = @(x, u) hdot_feas(x(1)*year, state(x(2:12)), constants.p1) ...
    + phi_feas(x(1)*year, state(x(2:12)), u, constants.p1)*dt ...
    + M2plus*dt + 1/2*M3plus*dt^2;
H_future1 = @(x, u) h_future1(x,u) + 1/(2*mu)*absSq(hdot_future1(x,u));

h_future2 = @(x, u) h_feas(x(1)*year, state(x(2:12)), constants.p2) ...
    + hdot_feas(x(1)*year, state(x(2:12)), constants.p2)*dt ...
    + 1/2*phi_feas(x(1)*year, state(x(2:12)), u, constants.p2)*dt^2 ...
    + 1/2*M2plus*dt^2 + 1/6*M3plus*dt^3;
hdot_future2 = @(x, u) hdot_feas(x(1)*year, state(x(2:12)), constants.p2) ...
    + phi_feas(x(1)*year, state(x(2:12)), u, constants.p2)*dt ...
    + M2plus*dt + 1/2*M3plus*dt^2;
H_future2 = @(x, u) h_future2(x,u) + 1/(2*mu)*absSq(hdot_future2(x,u));

func1 = @(x) h_omega_future(x, u_verify(x(1)*year, x(2:12), constants)); % there's no delta here because this is a first order constraint
func2 = @(x) (h_future1(x, u_verify(x(1)*year, x(2:12), constants)) + delta_2);
func3 = @(x) (H_future1(x, u_verify(x(1)*year, x(2:12), constants)) + Delta_2);
func4 = @(x) (h_future2(x, u_verify(x(1)*year, x(2:12), constants)) + delta_2);
func5 = @(x) (H_future2(x, u_verify(x(1)*year, x(2:12), constants)) + Delta_2);
func = @(x) -max([func1(x), func2(x), func3(x), func4(x), func5(x)])*1e4;
N = 100;
min_val = Inf;
min_res = [];

for i=1:min(N, N_max)
guess = randn(12,1).*upper;
try
[res, fval] = fmincon(func, guess, [], [], [], [], lower, upper, nonlcon);
if sum(nonlcon(res) > 0) == 0
    min_val = min(min_val, fval); if min_val == fval, min_res = res; end
end
catch e
disp(e)
end
end
min_val = min_val/1e4
if ~isempty(min_res), nonlcon(min_res), end

% If min_val is positive or zero for enough tests, then all the constraints together are
% jointly feasible.
% This section of the code can be troublesome at times and takes a long time to run, so 
% there is a degree of assumption that the number of cases tested so far is representative
% of the entire safe set. 

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
p1 = [1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];
p2 = [-sqrt(3)/2; 1/2; 0];
save('parameters.mat', 'Je', 'wheel_axis', 'wheel_limit', 'wheel_rate_limit',...
    'w_max', 'E_max', 'ctheta', 'dt', 'p1', 'p2', 'Z', 'torque', 'Jtot', 'Jw', ...
    'omega_Earth', 'O_Earth_from_Sol', 'delta_2', 'Delta_2', 'mu', ...
    'M2plus', 'M2minus', 'M3plus', 'M3minus', 'M1omega', 'M2omega', 'M2omega_lin');

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