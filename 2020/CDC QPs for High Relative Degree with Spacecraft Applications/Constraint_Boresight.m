function [A, b] = Constraint_Boresight(x, con, t, type)
% Returns the A, b, Aeq, and beq matrices for ensuring that the boresight requirement is
% achieved.
persistent mu ftilde_mu
if isempty(mu)
    file = load('InData/SimData.mat');
    mu = file.mu;
    ftilde_mu = file.ftilde_mu;
end
N = 0.1;
z_g = con.vector;
beta = con.scalar;
epsilon = con.epsilon;
t1 = con.t1;
t2 = con.t2;

rvec = x(1:3);
vvec = x(4:6);
yvec = z_g - rvec;
v = norm(vvec);
y = norm(yvec);
yhat = yvec/y;

mrp = x(7:9);
omega = x(10:12);
q = PtoQ(mrp);
e3 = [0; 0; 1];
what = RotQ(e3, q);
whatdot = cross(omega, what);

A1 = -1/y*what'*skew(yhat)^2;
A2 = yhat'*skew(what);

h = cos(beta) - what'*yhat;
hdot = -yhat'*whatdot - 1/y*what'*skew(yhat)^2*vvec;

if isequal(type, 'future')
    rhs = Quadratic_Future(t, t1, h, hdot);
elseif isequal(type, 'always')
    if h > 0
        rhs = -N;
    else
        rhs = Quadratic_Always(t, t2, h, hdot);
    end
else
    error('Temporal Logic was neither ''future'' or ''always''. Other cases not yet coded');
end
lambda = -norm(cross(yhat, cross(yhat, what)))*ftilde_mu/y;
b = rhs + 2/y*vvec'*skew(yhat)^2*whatdot + what'*( ...
    - 2*(yhat'*vvec)/y^2*vvec + 1/y^2*(3*(yhat'*vvec)^2 - v^2)*yhat) ...
    + yhat'*cross(omega, cross(omega, what)) + lambda;

A = [A1, A2];

if isequal(type, 'always')
%     tol = deg2rad(2);
%     if beta - acos(what'*yhat) > tol
    if h < -epsilon
        A = zeros(size(A));
        b = zeros(size(b));
    end
end
end