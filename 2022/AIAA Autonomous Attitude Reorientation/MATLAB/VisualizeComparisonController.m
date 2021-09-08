% This accessroy file is used to compare our quaternion notation against that in other
% papers, in particular Lee/Mesbahi 2010. The conclusion is that what we call the attitude
% quaternion is the conjugate of what this comparison paper calls the attitude quaternion.

global O_Earth_from_Sol omega_Earth
t = 0;
x1 = O_Earth_from_Sol*[cos(omega_Earth*t); sin(omega_Earth*t); 0]; % inertial frame
y1 = constants.p1;
ctheta1 = constants.ctheta;

M1 = [dot(x1, y1) - ctheta1, cross(x1, y1).';
    cross(x1, y1), x1*y1.' + y1*x1.' - (dot(x1,y1)+ctheta1)*eye(3)];

figure(1); clf;
PlotCone(get_s(t(1)), acos(constants.ctheta));
hold on;
[s1, s2, s3] = sphere(20);
surf(s1, s2, s3, 'FaceAlpha', 0);
xlabel 'x'; ylabel 'y'; zlabel 'z'; hold on;
axis equal;

for i=1:100
    q = randn(4,1);
    if q'*M1*q <= 0
        color = [0;1;0];
    else
        color = [1;0;0];
    end
    % My notation: QtoR(q) = R_{B/A} = O_{A/B} so QtoR(q)*p1_{B} = p1_{A}
%     p = QtoR(q)*constants.p1;
    
    % Their notation: though not thoroughly explained, it is clear from this experiment
    % that this is what needs to be done for the red dots to fall within the red cone.
    p = QtoR(q)'*constants.p1;
    plot3(p(1), p(2), p(3), 'o', 'Color', color, 'MarkerFaceColor', color);
end

% Conclusion: Everything in Lee/Mesbahi's paper is in terms of the conjugate of what we
% are calling the attitude quaternion.
% I have adjusted the matrices M1 and M2 in the ComparisonController.m file to convert my
% notation into their notation.