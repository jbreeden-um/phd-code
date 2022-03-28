function out = propagate_orbit(x,dt)
mu = 398600;
orbit = CartesianToKepler(x(1:3), x(4:6), mu);
M1 = nu_to_M(orbit.nu, orbit.e);
n = sqrt(mu/orbit.a^3);
M2 = M1 + n*dt;
nu2 = M_to_nu(M2, orbit.e);
out = KeplerToCartesian(orbit.a, orbit.e, orbit.i, orbit.Omega, orbit.omega, nu2, mu);
end

function E = nu_to_E(nu, e)
    E = atan2(sqrt(1-e^2)*sin(nu), e + cos(nu));
end

function M = nu_to_M(nu, e)
    E = nu_to_E(nu, e);
    M = E - e*sin(E);
end

function nu = M_to_nu(M, e)
    E = fzero(@(E) E - e*sin(E) - M, M);
    nu = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end

