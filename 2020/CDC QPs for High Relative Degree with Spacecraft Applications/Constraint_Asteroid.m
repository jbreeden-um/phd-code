function [A, b] = Constraint_Asteroid(x, con, t, type)
% Returns the A and b constraint matrices for ensuring avoidance is achieved.
% Results should be fed into QP as A \leq b
persistent mu ftilde_mu Eros
if isempty(mu)
    file = load('InData/SimData.mat');
    mu = file.mu;
    ftilde_mu = file.ftilde_mu;
    file = load('InData/Eros_Shape.mat');
    Eros = (file.vertices)'; % Note that Eros's shape model file is in units of km
end
N = 10;
rhobar = con.scalar;
epsilon = con.epsilon;
t1 = con.t1;
t2 = con.t2;

rvec = x(1:3);
vvec = x(4:6);
v = norm(vvec);

vecs = rvec - Eros;
dists = vecnorm(vecs);
[~,index] = min(dists);
z_s = Eros(:,index);

h = rhobar^2 - (rvec-z_s)'*(rvec-z_s);
hdot = 2*(z_s-rvec)'*vvec;

A1 = -2*(rvec-z_s)';
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
lambda = -2*norm(rvec - z_s)*ftilde_mu;
b = rhs + 2*v^2 + lambda;

A = [A1, zeros(1,3)]; % We don't care about the attitude input here
if isequal(type, 'always')
%     tol = 2;
%     if norm(rvec-z_s) > rhobar+tol
    if h < -epsilon
        A = zeros(size(A));
        b = zeros(size(b));
    end
end
end