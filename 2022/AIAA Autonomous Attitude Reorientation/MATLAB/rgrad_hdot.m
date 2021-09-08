function out = rgrad_hdot(t,x)
global O_Earth_from_Sol omega_Earth

s = O_Earth_from_Sol*[cos(omega_Earth*t); sin(omega_Earth*t); 0];
sdot = O_Earth_from_Sol*[-omega_Earth*sin(omega_Earth*t); omega_Earth*cos(omega_Earth*t); 0];
r = x(1:3);
w = x(4:6);

out = [sdot' + s'*skew(w), -s'*skew(r), zeros(1,4)];
end