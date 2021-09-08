function out = phi(t,x,u,p)
global O_Earth_from_Sol omega_Earth Z Jtot wheel_axis Jw
s = QtoR(x(1:4))'*O_Earth_from_Sol*[cos(omega_Earth*t); sin(omega_Earth*t); 0];
sdot = QtoR(x(1:4))'*O_Earth_from_Sol*[-omega_Earth*sin(omega_Earth*t); omega_Earth*cos(omega_Earth*t); 0];
sddot = QtoR(x(1:4))'*O_Earth_from_Sol*[-omega_Earth^2*cos(omega_Earth*t); -omega_Earth^2*sin(omega_Earth*t); 0];
r = p;
w = x(5:7);
W = x(8:11);
wdot = Z(1:3,:)*[-cross(w,Jtot*w + wheel_axis*Jw*W); u];

out = sddot'*r + 2*sdot'*skew(w)*r + s'*skew(w)^2*r - s'*skew(r)*wdot;
end