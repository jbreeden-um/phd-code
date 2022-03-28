function out = propagate_debris(t)
% Assume circular orbit
global debris_ke
mu = 398600;
n = sqrt(mu/debris_ke.a^3);
nu2 = debris_ke.nu0 + n*t;
out = KeplerToCartesian(debris_ke.a, debris_ke.e, debris_ke.i, debris_ke.O, debris_ke.o, nu2, mu);
end