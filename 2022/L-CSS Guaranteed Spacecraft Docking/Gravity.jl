module Gravity

using LinearAlgebra
include("Utilities.jl");

#################################################

r_Ceres = 476000;
mu_Ceres = 6.26325e10;

function gravity_Ceres(t, r)
    return -mu_Ceres*r/norm(r)^3;
end

function partial_gravity_Ceres(t, r, vc=[0;0;0])
    return 3*mu_Ceres*r*r'/norm(r)^5*(-vc) - mu_Ceres/norm(r)^3*(-vc);
end

I3 = [1.0 0.0 0.0;
      0.0 1.0 0.0;
      0.0 0.0 1.0];
function grad_gravity_Ceres(t, r)
    return 3*mu_Ceres*r*r'/norm(r)^5 - mu_Ceres*I3/norm(r)^3;
end

function obstacle_Ceres(t, args)
    rc = [0;0;0];
    vc = [0;0;0];
    uc = [0;0;0];
    udotc = [0;0;0];
    return rc, vc, uc, udotc;
end

#################################################

mu_Earth = 398600e9;
R_Earth = 6378e3;
altitude_Earth = 400e3;
n = sqrt(mu_Earth / (R_Earth + altitude_Earth)^3);

function HCW(t, x)
    return [0 0 1 0; 0 0 0 1; 3*n^2 0 0 2*n; 0 0 -2*n 0]*x;
end

end