function u = ComparisonControllerSMC(t, x)
global s_target
persistent constants k2 k3 alpha_V Jinv mapping
if t<0
    constants = load('parameters.mat');
    
    % Gains
    k2 = 3e4*max(eig(constants.Jtot));
    k3 = 0.1*max(eig(constants.Jtot));
    alpha_V = 1/50;
    Jinv = inv(constants.Jtot);
    mapping = pinv(constants.Z(1:3,4:7));
end
global outdata

q = x(1:4);
w = x(5:7);

% Desired Quaternion Calculation
sB = QtoR(x(1:4))'*s_target; % command vector in body frame
rB = constants.p1; % current pointing vector in body frame
orth = cross(rB, sB); % vector to rotate about in body frame
orth = orth/norm(orth);
alpha = dot(sB, rB);
if abs(alpha) > 1
    alpha = 1*sign(alpha);
end
% Note that if one limits alpha to 0.02 as we did for the linear
% controller, then this approach will be poorly posed.
alpha = acos(alpha); % rotation angle
dq = [cos(alpha/2); sin(alpha/2)*orth]; % quaternion from body frame to command frame
qr = QxQ(x(1:4), dq);
qe = QxQ(QConj(qr), q);

% Set up exlcusion zones
x1 = get_s(t); % inertial frame
y1 = constants.p1;
y2 = constants.p2;

global ctheta
ctheta1 = ctheta;
ctheta2 = ctheta;

% Gradients
% grad_V = 2*(q-qr)*alpha_V; % Use this to ignore the constraints
grad_V = ComparisonDerivativesSMC(q, qr, x1, y1, y2, ctheta1, ctheta2, alpha_V); % Use this to enforce the constraints

k = 0.01; % This is the largest that one can make k without things breaking
s = w + k*qe(2:4);
s_max = constants.w_max - k; % Note that w max refers to the largest axis (slightly different from the energy based constraint in the paper)
if sum(abs(s) > s_max)
    disp('s is too large');
end
Psi = diag(s_max^2 - s.^2);
f = -skew(w)*constants.Jtot*w + k/2*(skew(qe(2:4)) + qe(1)*eye(3))*w;

term1 = QxQ(QConj(grad_V), q);
Gamma = Psi*Jinv;
u0 = -Gamma*(k2*s - k3*term1(2:4)) - f;
    
u = mapping*u0;

if norm(u) > 0.1
    disp('Look at this');
end

u = limit(u, -constants.wheel_limit*[1;1;1;1], constants.wheel_limit*[1;1;1;1]);

outdata.u = u;

if sum(isnan(u))
    disp 'NaNs';
end

end