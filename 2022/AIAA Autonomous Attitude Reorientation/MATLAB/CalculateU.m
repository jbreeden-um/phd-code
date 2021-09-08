function u = CalculateU(t, x)
global s_target use_lyapunov use_cbf
persistent constants
if t==0
    constants = load('parameters.mat');
end
global outdata

% CBF Conditions
[H0, A0, b0] = ConstraintE(t, x, constants, 3);
[A1, b1] = ConstraintQ(t, x, constants, constants.p1, 1);
[A2, b2] = ConstraintQ(t, x, constants, constants.p2, 2);

nonlcon_omega = @(u) u'*H0*u + A0*u - b0;

scale = 1e5;
lower = -constants.wheel_limit*[1;1;1;1]*scale;
upper =  constants.wheel_limit*[1;1;1;1]*scale;

sB = QtoR(x(1:4))'*s_target; % command vector in body frame
rB = constants.p1; % current pointing vector in body frame
orth = cross(rB, sB); % vector to rotate about in body frame
orth = orth/norm(orth);
alpha = dot(sB, rB);
if abs(alpha) > 1
    alpha = 1*sign(alpha);
end
alpha = min(acos(alpha), 0.2); % rotation angle
dq = [cos(alpha/2); sin(alpha/2)*orth]; % quaternion from body frame to command frame
omega = x(5:7);
domega = [0;0;0] - omega;

kp = 0.1;
kd = 0.5;
u0 = kp*dq(2:4) + kd*domega;
% This controller has no integral term, so it will not work for fine
% pointing control. It is intended for reorientation.
if ~use_lyapunov
    wheels0 = pinv(constants.Z(1:3,4:7))*u0;
else
    wheels0 = ComparisonController(t,x);
end
if ~use_cbf
    u = wheels0;
    outdata.u = u;
    return;
end

F = -2*wheels0;
A = [A1; A2];
b = [b1;b2];

try
    [u, fval, status] = fmincon(@(u) (u'*u/scale^2 + F'*u/scale)*1e10, ...
        [0;0;0;0], A, b*scale, [], [], lower, upper, ...
        @(x) assemble_constraint(x, @(y) nonlcon_omega(y/scale)*1e6), ...
        optimset('Display', 'off', 'TolX', 1e-8));
    u = u/scale;
    outdata.u = u;
catch e
    disp(e);
end

if status < 1
    disp(['Something went wrong. Status: ' num2str(status)]);
end

if sum(isnan(u))
    disp 'NaNs';
end

end

function out = limit(x, min, max)
out = x;
for i=1:length(x)
    if out(i) < min(i)
        out(i) = min(i);
    elseif out(i) > max(i)
        out(i) = max(i);
    end
end
end