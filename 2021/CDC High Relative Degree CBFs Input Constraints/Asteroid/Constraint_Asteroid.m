function [A, b] = Constraint_Asteroid(x)
% Returns the A and b constraint matrices for ensuring avoidance is achieved.
% Results should be fed into QP as A \leq b
persistent a_max mu ftilde
if isempty(a_max)
    file = load('InData/SimData.mat');
    u_max = file.u_max;
    file = load('InData/Eros_GravApprox.mat');
    mu = file.mu;
    ftilde = file.ftilde;
    a_max = u_max - file.fhatmax;
end

rvec = x(1:3);
vvec = x(4:6);
rho = 1;

[A, b] = FindOtherPoints(rvec, vvec, rho, a_max, mu, ftilde);
end

function [A, b] = FindOtherPoints(rvec, vvec, rho, a_max, mu, ftilde)
% Calculates a subset of active constraints and computes the A and b matrices for each
% one.
global is_error outdata plates_to_consider
persistent points normals
if isempty(points)
    file = load('InData/Eros_Points.mat');
    points = file.points;
    normals = file.normals;
end

v = norm(vvec);
if v > 0
    vhat = vvec(:)/v;
else
    vhat = vvec(:);
end

diffs = points - rvec(:)';
parallel = (diffs*vhat)*vhat';
distances = vecnorm(parallel,2,2);
projs = diffs - parallel;
margin = sind(10);
to_consider = dot(normals, diffs, 2) <= 0 & vecnorm(projs,2,2) <= rho + distances*margin & diffs*vvec >= 0;
plates_to_consider = to_consider;

n = sum(to_consider);
rc_all = points(to_consider, :)';
outdata.rc_all = rc_all';

drvec = -diffs(to_consider,:);
dvvec = vvec;
r = vecnorm(drvec,2,2);
v = norm(dvvec);

h = rho - r;
hdot = -(drvec*dvvec)./r;
fhat_mu = -mu/norm(rvec)^3*rvec(:);

LgLf = -drvec./r;
Lf2 = ((drvec*dvvec).^2 - r.^2*v^2)./r.^3 - drvec*fhat_mu./r;
if sum(hdot < 0) > 0, disp 'Check your math'; end
H = h + hdot.^2/(a_max*2);
if sum(H > 0) > 0, is_error = 1; end

lambda = hdot*ftilde/a_max; % gravity error compensation term
alpha = -H; % straight line of slope 1
A = hdot.*LgLf/a_max;
b = alpha - hdot.*(1+Lf2/a_max) - lambda;

outdata.n_avoid = n;
[~,i] = min(vecnorm(diffs,2,2));
outdata.drc_alg = diffs(i,:)';
outdata.H_max = max(H);
outdata.h_max = max(h);
end