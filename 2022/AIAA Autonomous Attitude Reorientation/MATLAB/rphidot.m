function out = rphidot(t,x,u,v)
global O_Earth_from_Sol omega_Earth Z Jtot wheel_axis Jw
s = O_Earth_from_Sol*[cos(omega_Earth*t); sin(omega_Earth*t); 0];
sdot = O_Earth_from_Sol*[-omega_Earth*sin(omega_Earth*t); omega_Earth*cos(omega_Earth*t); 0];
sddot = O_Earth_from_Sol*[-omega_Earth^2*cos(omega_Earth*t); -omega_Earth^2*sin(omega_Earth*t); 0];
sdddot = O_Earth_from_Sol*[omega_Earth^3*sin(omega_Earth*t); -omega_Earth^3*cos(omega_Earth*t); 0];
r = x(1:3);
w = x(4:6);
W = x(7:10);
wdot = Z(1:3,:)*[-cross(w,Jtot*w + wheel_axis*Jw*W); u];
Wdot_v = Z(4:7,:)*[-cross(w,Jtot*w + wheel_axis*Jw*W) + v; u];

wdot_v = Z(1:3,:)*[-cross(w,Jtot*w + wheel_axis*Jw*W) + v; u];
wddot_v = Z(1:3,1:3)*(-cross(wdot_v,Jtot*w + wheel_axis*Jw*W) - cross(w,Jtot*wdot_v + wheel_axis*Jw*Wdot_v));

out = sdddot'*r + 3*sddot'*skew(w)*r + 3*sdot'*skew(w)^2*r - sdot'*skew(r)*wdot - 2*sdot'*skew(r)*wdot_v ...
    + s'*skew(wdot)*skew(w)*r + s'*skew(wdot_v)*skew(w)*r + s'*skew(w)*skew(wdot_v)*r + s'*skew(w)^3*r ...
    - s'*skew(r)*wddot_v;
end