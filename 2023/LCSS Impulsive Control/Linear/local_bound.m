function mag = local_bound(x)
% Given a current x value, this function returns the maximum value of norm(x) that can
% occur anywhere in the current sublevel set of the Lyapunov function. That is, if the
% Lyapunov function is nonincreasing, we can use this to construct psi_h and psi_v
global P
level = x'*P*x;
r = eig(P);
e_min = min(r);
mag = sqrt(level/e_min);
end