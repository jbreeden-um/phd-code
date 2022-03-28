function out = lane2(z2)
global sim_case
if sim_case==1
    out = lane2_parallel(z2);
elseif sim_case==2
    out = lane2_left(z2);
else
    error('Unknown simulation case.');
end
end

function out = lane2_parallel(z2)
out = [0; z2] + [1.5; 0];
end

function out = lane2_left(z2)
arc_length = pi/2*4.5;
if z2 <= -3
    out = [0; z2] + [1.5; 0];
elseif z2 <= -3+arc_length
    theta = (z2 + 3)/arc_length*pi/2;
    out = [-3 + 4.5*cos(theta); -3 + 4.5*sin(theta)];
else
    out = [-3 - (z2+3-arc_length); 1.5];
end
end