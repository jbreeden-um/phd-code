function u = ComparisonControllerMPC(t, x)
global s_target outdata ctheta
persistent constants mapping mpc alpha qr
if t<0
    constants = load('parameters.mat');
    mapping = pinv(constants.Z(1:3,4:7));

    % Parameters
    mpc = struct();
    mpc.h = constants.dt;
    mpc.alpha = constants.wheel_limit*min(svd(constants.wheel_axis));
    mpc.beta1 = -constants.ctheta; % inclusion angle
    mpc.v1 = constants.p1; % body vector
    mpc.w1 = -get_s(t); % inertial vector (negative of the vector to be excluded)

    mpc.beta2 = -ctheta;
    mpc.v2 = constants.p2;
    mpc.w2 = -get_s(t);

    mpc.J = constants.Jtot;
    mpc.Jd = 0.5*trace(mpc.J)*eye(3)-mpc.J;

    mpc.Q1 = 1e-2*eye(3);
    mpc.Q2 = 3.8e1*eye(3);
    mpc.Q3 = 1e2*eye(3);

    mpc.P1 = 1e-2*eye(3);
    mpc.P2 = 3.8e1*eye(3);
    
%     Barrier = load('Results/Barrier.mat');
%     mpc.R_target = QtoR(Barrier.x(end,1:4)');
    alpha = 1;
end

% Desired Quaternion Calculation
if alpha > 0.1
    sB = QtoR(x(1:4))'*s_target; % command vector in body frame
    rB = constants.p1; % current pointing vector in body frame
    orth = cross(rB, sB); % vector to rotate about in body frame
    orth = orth/norm(orth);
    alpha = dot(sB, rB);
    if abs(alpha) > 1
        alpha = 1*sign(alpha);
    end
    alpha = acos(alpha);
    dq = [cos(alpha/2); sin(alpha/2)*orth]; % quaternion from body frame to command frame
    % If we do not eventually fix the target, this particular MPC
    % implementation will never converge
    qr = QxQ(x(1:4), dq);
end

% Notations for this file:
% R = O_{inertial/body}
% h = cos(beta) - dot((R*v), w) \leq 0 

mpc.w1 = -get_s(t);
mpc.w2 = -get_s(t);
mpc.R_target = QtoR(qr);

R = QtoR(x(1:4));
AM = constants.Jtot*x(5:7);

n = 5;
Tp = n*mpc.h;
mu = [1e10;1e3];
par = [0;0;0;0;0;0];
eta0 = [0;0;0];
dAM0 = [0;0;0];
    
paropt = mpc_solver(R,AM,Tp,mu,[0;0],eta0,dAM0,par,mpc);
    
[u0,new_R,new_AM,F,M,O,p1,p2] = mpc_func3(R,AM,mpc.h,mu,paropt,mpc);
    
u = mapping*u0;

u = limit(u, -constants.wheel_limit*[1;1;1;1], constants.wheel_limit*[1;1;1;1]);

outdata.u = u;

if sum(isnan(u))
    disp 'NaNs';
end

end