function out = GradH_func(x)
global constants
sigma = constants.sigma;
r = x(1:2); phi = x(3);
out = [-((2*r(1)-(2*sigma*wrap_angle(-atan2(r(2),r(1))+phi)*r(2))/(r(1)^2+r(2)^2))/(2*sqrt(-sigma ...
        *(-(pi^2/4)+wrap_angle(-atan2(r(2),r(1))+phi)^2)+r(1)^2+r(2)^2))), ...
    -((2*r(2)+(2*sigma*wrap_angle(-atan2(r(2),r(1))+phi)*r(1))/(r(1)^2+r(2)^2))/(2*sqrt(-sigma ...
        *(-(pi^2/4)+wrap_angle(-atan2(r(2),r(1))+phi)^2)+r(1)^2+r(2)^2))), ...
    (sigma*wrap_angle(-atan2(r(2),r(1))+phi))/sqrt(-sigma*(-(pi^2/4)+wrap_angle(-atan2(r(2),r(1))+phi)^2)+r(1)^2+r(2)^2)];
end