function [r, v] = obstacle_location(index,t)
persistent oe1 oe2
mu = 398600e9;
Re = 6378e3;
alt = 600e3;
n = sqrt(mu/(Re+alt)^3);
R = [cos(n*t), -sin(n*t); sin(n*t), cos(n*t)]; R = blkdiag(R, R);
if index==1
    if t==0
        d0 = R*[0; -7.5e3; 0; 0];
        x0 = d0 + get_center(t) + [0; 0; [eye(2), zeros(2,1)]*cross([0;0;n], [d0(1:2); 0])];
        oe1 = CartesianToKepler([x0(1:2); 0], [x0(3:4); 0], mu);
    end
    oe = oe1;
elseif index==2
    if t==0
        d0 = R*[0; -2.5e3; 0; 0];
        x0 = d0 + get_center(t) + [0; 0; [eye(2), zeros(2,1)]*cross([0;0;n], [d0(1:2); 0])];
        oe2 = CartesianToKepler([x0(1:2); 0], [x0(3:4); 0], mu);
    end
    oe = oe2;
elseif index==3
    d0 = R*[6e3; -5e3; 0; 0];
    x0 = d0 + get_center(t) + [0; 0; [eye(2), zeros(2,1)]*cross([0;0;n], [d0(1:2); 0])];
elseif index==4
    d0 = R*[-6e3; -5e3; 0; 0];
    x0 = d0 + get_center(t) + [0; 0; [eye(2), zeros(2,1)]*cross([0;0;n], [d0(1:2); 0])];
else
    error('Unreognized obstacle location');
end
if index==1 || index==2
    n = sqrt(mu/oe.a^3);
    x = KeplerToCartesian(oe.a, oe.e, oe.i, oe.Omega, oe.omega, oe.nu + n*t, mu);
    r = x(1:2);
    v = x(4:5);
else
    % Note that because h_derivations assumes that all the obstacles are follow real
    % \ddot{r}=-\mu*r/norm(r)^3 dynamics, the results of this simulation are not strictly
    % correct with regard to obstacles 3 and 4. However, otherwise, these obstacles
    % quickly leave the sphere of influence.
    r = x0(1:2);
    v = x0(3:4);
end
end