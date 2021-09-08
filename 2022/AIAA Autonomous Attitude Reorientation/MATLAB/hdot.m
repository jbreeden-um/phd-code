function out = hdot(t,x,p)
global O_Earth_from_Sol omega_Earth

s = QtoR(x(1:4))'*O_Earth_from_Sol*[cos(omega_Earth*t); sin(omega_Earth*t); 0];
sdot = QtoR(x(1:4))'*O_Earth_from_Sol*[-omega_Earth*sin(omega_Earth*t); omega_Earth*cos(omega_Earth*t); 0];
r = p;
w = x(5:7);

out = sdot'*r + s'*skew(w)*r;
end
