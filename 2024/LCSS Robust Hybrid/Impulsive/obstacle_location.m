function [y, oe] = obstacle_location(index,t)
persistent oe1 oe2 oe3 oe4 oe5 oe6
mu = 398600e9;
Re = 6378e3;
alt = 600e3;
n = sqrt(mu/(Re+alt)^3);
[xc, theta, oe_c] = get_center(t);
R = [cos(theta), -sin(theta); sin(theta), cos(theta)]; R = blkdiag(R, R);
if index==1
    if t==0
        d0 = R*[0; -10e3; 0; 0];
        x0 = d0 + xc + [0; 0; [eye(2), zeros(2,1)]*cross([0;0;n], [d0(1:2); 0])];
        oe1 = CartesianToKepler([x0(1:2); 0], [x0(3:4); 0], mu);
        oe1.a = oe_c.a;
    end
    oe = oe1;
elseif index==2
    if t==0
        d0 = R*[0; -4e3; 0; 0];
        x0 = d0 + xc + [0; 0; [eye(2), zeros(2,1)]*cross([0;0;n], [d0(1:2); 0])];
        oe2 = CartesianToKepler([x0(1:2); 0], [x0(3:4); 0], mu);
        oe2.a = oe_c.a;
    end
    oe = oe2;
elseif index==3
    if t==0
        d0 = R*[3e3; -5e3; 0; 0];
        x0 = d0 + xc + [0; 0; [eye(2), zeros(2,1)]*cross([0;0;n], [d0(1:2); 0])];
        oe3 = CartesianToKepler([x0(1:2); 0], [x0(3:4); 0], mu);
        oe3.a = oe_c.a;
    end
    oe = oe3;
elseif index==4
    if t==0
        d0 = R*[-3e3; -5e3; 0; 0];
        x0 = d0 + xc + [0; 0; [eye(2), zeros(2,1)]*cross([0;0;n], [d0(1:2); 0])];
        oe4 = CartesianToKepler([x0(1:2); 0], [x0(3:4); 0], mu);
        oe4.a = oe_c.a;
    end
    oe = oe4;
elseif index==5
    if t==0
        delta_t = 1600;
        [xc, theta] = get_center(delta_t);
        R = [cos(theta), -sin(theta); sin(theta), cos(theta)]; R = blkdiag(R, R);
        d0 = R*[1.3e3; -5e3; 0; 0];
        x0 = d0 + xc + [0; 0; [eye(2), zeros(2,1)]*cross([0;0;n], [d0(1:2); 0])];
        oe5 = CartesianToKepler([x0(1:2); 0], [x0(3:4); 0], mu);
        oe5.a = oe_c.a;
        E = 2*atan(sqrt((1-oe5.e)/(1+oe5.e))*tan(oe5.nu/2));
        M = E - oe5.e*sin(E);
        n = sqrt(mu/oe5.a^3); 
            % This n is different from the n computed at the top of the file. This is okay
            % since the obstacles are all artificial anyways.
        delta_M = -n*delta_t;
        M0 = M + delta_M;
        oe5.nu = convert_M_to_nu(M0,oe5.e);
    end
    oe = oe5;
elseif index==6
    if t==0
        delta_t = 1600;
        [xc, theta] = get_center(delta_t);
        R = [cos(theta), -sin(theta); sin(theta), cos(theta)]; R = blkdiag(R, R);
        d0 = R*[-1.7e3; -4e3; 0; 0];
        x0 = d0 + xc + [0; 0; [eye(2), zeros(2,1)]*cross([0;0;n], [d0(1:2); 0])];
        oe6 = CartesianToKepler([x0(1:2); 0], [x0(3:4); 0], mu);
        oe6.a = oe_c.a;
        E = 2*atan(sqrt((1-oe6.e)/(1+oe6.e))*tan(oe6.nu/2));
        M = E - oe6.e*sin(E);
        n = sqrt(mu/oe6.a^3);
        delta_M = -n*delta_t;
        M0 = M + delta_M;
        oe6.nu = convert_M_to_nu(M0,oe6.e);
    end
    oe = oe6;
else
    error('Unreognized obstacle location');
end
n = sqrt(mu/oe_c.a^3);
E0 = 2*atan(sqrt((1-oe.e)/(1+oe.e))*tan(oe.nu/2));
M0 = E0 - oe.e*sin(E0);
M = M0 + n*t;
nu = convert_M_to_nu(M, oe.e);
x = KeplerToCartesian(oe.a, oe.e, oe.i, oe.Omega, oe.omega, nu, mu);
y = x([1;2;4;5]);
end