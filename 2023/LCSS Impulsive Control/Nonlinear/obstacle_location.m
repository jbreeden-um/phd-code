function [r, v] = obstacle_location(index,t)
global A_sys use_moving_obstacles
if index==1
    x0 = [0; -7.5e3; 0; 0];
elseif index==2
    x0 = [0; -2.5e3; 0; 0];
elseif index==3
    x0 = [6e3; -5e3; 1; 0];
elseif index==4
    x0 = [-6e3; -5e3; 1; 0];
else
    error('Unreognized obstacle location');
end
Phi = expm(A_sys*t*use_moving_obstacles);
r = Phi(1:2,:)*x0;
v = Phi(3:4,:)*x0;
end