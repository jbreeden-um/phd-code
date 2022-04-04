function varargout = CartesianToKepler(rvec, vvec, mu)
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
if nhat(2)<0, Omega = 2*pi-Omega; end;
omega = acos(dot(nhat, ehat));
if ehat(3)<0, omega = 2*pi-omega; end;
nu = acos(dot(ehat, rvec/r));
if dot(rvec, vvec)<0, nu = 2*pi-nu; end;
if nargout==6
    varargout = {a, e, i, Omega, omega, nu};
elseif nargout==1 || nargout==0
    varargout = {struct('a',a,'e',e,'i',i,'Omega',Omega,'omega',omega,'nu',nu)}; 
else
    disp 'Provide either 1 or 6 output arguments';
end