function varargout = KeplerToCartesian(a, e, i, Omega, omega, nu, mu)
r = a*(1-e^2)/(1+e*cos(nu));
rvecP = [r*cos(nu); r*sin(nu); 0];
vvecP = sqrt(mu/(a*(1-e^2)))*[-sin(nu); (e+cos(nu)); 0];
RP2I = [cos(Omega), -sin(Omega), 0; sin(Omega), cos(Omega), 0; 0, 0, 1]...
    *[1, 0, 0; 0, cos(i), -sin(i); 0, sin(i), cos(i)]...
    *[cos(omega), -sin(omega), 0; sin(omega), cos(omega), 0; 0, 0, 1];
rvecI = RP2I*rvecP;
vvecI = RP2I*vvecP;
if nargout==2
    varargout = {rvecI, vvecI};
elseif nargout==1 || nargout==0
    varargout = {[rvecI; vvecI]};
end
end