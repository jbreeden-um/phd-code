function u = ComparisonController(t, x)
global s_target O_Earth_from_Sol omega_Earth k_compare
persistent constants
if t==0
    constants = load('parameters.mat');
end
global outdata

% Desired Quaternion Calculation
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
qr = QxQ(x(1:4), dq);

% Set up exlcusion zones
x1 = O_Earth_from_Sol*[cos(omega_Earth*t); sin(omega_Earth*t); 0]; % inertial frame
y1 = constants.p1;
y2 = constants.p2;

global ctheta
ctheta1 = ctheta;
ctheta2 = ctheta;

% Nonlinear Backstepping Parameters
alpha = 0.75;
beta = 8;
k = k_compare;

% Gradients
q = x(1:4);
[grad_Vq, hess_Vq] = ComparisonDerivatives(q, qr, x1, y1, y2, ctheta1, ctheta2, k);
    % Note: the above function essentially takes the conjugate of q in the definition of
    % the avoidance criteria. This is due to a notational mismatch which only affects the
    % avoidance criteria and not the convergence criteria. This is handled entirely in the
    % above function and does not need to be considered in any of the remaining math.
wc = -2*QxQ(QConj(q), grad_Vq);
w = x(5:7);
z = alpha*atan(beta*(w - wc(2:4))); % the arctangent is applied elementwise

C2 = 1/(alpha*beta)*(eye(3) + beta^2*diag(w - wc(2:4))^2);
qdot = 1/2*QxQ(q, [0; w]); % expressed in body frame
wcdot = -2*(QxQ(QConj(qdot), grad_Vq) + QxQ(QConj(q), hess_Vq*qdot));
term3 = QxQ(QConj(grad_Vq), q);
u0 = constants.Jtot*wcdot(2:4) - skew(w)*constants.Jtot*w + 1/2*C2*term3(2:4) - z;
u = pinv(constants.Z(1:3,4:7))*u0;

u = limit(u, -constants.wheel_limit*[1;1;1;1], constants.wheel_limit*[1;1;1;1]);

outdata.u = u;

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