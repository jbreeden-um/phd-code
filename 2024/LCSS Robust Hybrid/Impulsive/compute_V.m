function V = compute_V(t,x1,get_P)
persistent P
if t==0
    mu = 398600e9;
    [~, ~, oe] = get_center(0);
    n = sqrt(mu/oe.a^3);
    A_sys = [0,     0, 1,    0;
             0,     0, 0,    1;
             3*n^2, 0, 0,    2*n;
             0,     0, -2*n, 0];
    B_sys = [zeros(2); eye(2)];
    R = eye(2)*100;
    Q = eye(4);
    [~, P, ~] = lqr(A_sys,B_sys,Q,R);
    if nargin==3 && get_P
        V = P;
        return
    end
end
xc = get_center(t);
V = (xc(1) - x1(1))*(P(1,1)*(xc(1) - x1(1)) + P(1,2)*(xc(2) - x1(2)) + P(1,3)*(xc(3) - x1(3)) + P(1,4)*(xc(4) - x1(4))) + (xc(2) - x1(2))*(P(1,2)*(xc(1) - x1(1)) + P(2,2)*(xc(2) - x1(2)) + P(2,3)*(xc(3) - x1(3)) + P(2,4)*(xc(4) - x1(4))) + (xc(3) - x1(3))*(P(1,3)*(xc(1) - x1(1)) + P(2,3)*(xc(2) - x1(2)) + P(3,3)*(xc(3) - x1(3)) + P(3,4)*(xc(4) - x1(4))) + (xc(4) - x1(4))*(P(1,4)*(xc(1) - x1(1)) + P(2,4)*(xc(2) - x1(2)) + P(3,4)*(xc(3) - x1(3)) + P(4,4)*(xc(4) - x1(4)));
end