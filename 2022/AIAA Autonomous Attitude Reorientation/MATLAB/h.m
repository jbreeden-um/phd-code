function out = h(t,x,p)
global O_Earth_from_Sol omega_Earth ctheta
r = p;
s = QtoR(x(1:4))'*O_Earth_from_Sol*[cos(omega_Earth*t); sin(omega_Earth*t); 0];

out = dot(s,r) - ctheta;
end