function out = get_s(t)
global O_Earth_from_Sol omega_Earth
out = O_Earth_from_Sol*[cos(omega_Earth*t); sin(omega_Earth*t); 0];
end
