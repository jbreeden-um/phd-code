function out = rhddot(t,x,u,v)
global O_Earth_from_Sol omega_Earth Z Jtot wheel_axis Jw
s = O_Earth_from_Sol*[cos(omega_Earth*t); sin(omega_Earth*t); 0];
sdot = O_Earth_from_Sol*[-omega_Earth*sin(omega_Earth*t); omega_Earth*cos(omega_Earth*t); 0];
sddot = O_Earth_from_Sol*[-omega_Earth^2*cos(omega_Earth*t); -omega_Earth^2*sin(omega_Earth*t); 0];
r = x(1:3);
w = x(4:6);
W = x(7:10);
wdot = Z(1:3,:)*[-cross(w,Jtot*w + wheel_axis*Jw*W)+v; u];

out = sddot'*r + 2*sdot'*skew(w)*r + s'*skew(w)^2*r - s'*skew(r)*wdot;
end