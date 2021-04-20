function [A, b] = Constraint_Proximity(x, con, t, type)
% Returns the A and b constraint matrices for ensuring proximity is achieved.
% Results should be fed into QP as A \leq b
persistent mu ftilde_mu
if isempty(mu)
    file = load('InData/SimData.mat');
    mu = file.mu;
    ftilde_mu = file.ftilde_mu;
end
N = 10;
z_g = con.vector;
rho = con.scalar;
epsilon = con.epsilon;
t1 = con.t1;
t2 = con.t2;

rvec = x(1:3);
vvec = x(4:6);
v = norm(vvec);

h = (z_g-rvec)'*(z_g-rvec) - rho^2;
hdot = -2*(z_g-rvec)'*vvec;

A1 = 2*(rvec-z_g)';
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
lambda = -2*norm(rvec - z_g)*ftilde_mu;
b = rhs - 2*v^2 + lambda;

A = [A1, zeros(1,3)]; % We don't care about the attitude input here

if isequal(type, 'always')
%     tol = 2;
%     if norm(z_g-rvec) < rho-tol
    if h < -epsilon
        A = zeros(size(A));
        b = zeros(size(b));
    end
end
end