
using LinearAlgebra

function absSq(x)
    return x*abs(x)
end

function KeplerToCartesian(a, e, i, Omega, omega, nu, mu)
    r = a*(1-e^2)/(1+e*cos(nu));
    rvecP = [r*cos(nu); r*sin(nu); 0];
    vvecP = sqrt(mu/(a*(1-e^2)))*[-sin(nu); (e+cos(nu)); 0];
    RP2I = [cos(Omega) -sin(Omega) 0; sin(Omega) cos(Omega) 0; 0 0 1] *
        [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)] *
        [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
    rvecI = RP2I*rvecP;
    vvecI = RP2I*vvecP;
    return rvecI, vvecI
end

function CartesianToKepler(rvec, vvec, mu)
    r = norm(rvec);
    v = norm(vvec);
    evec = ((v^2-mu/r)*rvec - dot(rvec,vvec)*vvec)/mu;
    e = norm(evec);
    ehat = evec/e;
    a = -mu/2*(v^2/2-mu/r)^(-1);
    hvec = cross(rvec, vvec);
    hhat = hvec/norm(hvec);
    i = acos(dot(hhat, [0;0;1]));
    nvec = cross([0;0;1], hvec);
    nhat = nvec/norm(nvec);
    Omega = acos(dot(nhat, [1;0;0]));
    if nhat[2]<0; Omega = 2*pi-Omega; end;
    omega = acos(dot(nhat, ehat));
    if ehat[3]<0; omega = 2*pi-omega; end;
    nu = acos(dot(ehat, rvec/r));
    if dot(rvec, vvec)<0; nu = 2*pi-nu; end;
    return a, e, i, Omega, omega, nu
end

"""
skew Returns skew symmetric matrix of v
* Returns the matrix by which one can premultiply an arbitrary vector to obtain the cross product of the input vector and the arbitary vector
"""
function skew(v)
    return [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0];
end