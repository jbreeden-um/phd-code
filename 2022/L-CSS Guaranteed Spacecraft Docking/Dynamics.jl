
include("Gravity.jl");

function UpdateX_Ceres!(xdot, x, p, t)
    Controller = p[1];
    Gravity = p[2];
    MatchedDisturbance = p[3];
    UnmatchedDisturbance = p[4];
    
    u = Controller(t,x);
    w_u = MatchedDisturbance(t);
    w_x = UnmatchedDisturbance(t);

    xdot[1:3] = x[4:6] + w_x;
    xdot[4:6] = Gravity(t, x[1:3]) + w_u + u;
end

function UpdateX_HCW!(xdot, x, p, t)
    Controller = p[1];
    Gravity = p[2];
    MatchedDisturbance = p[3];
    UnmatchedDisturbance = p[4];
    
    u = Controller(t,x);
    w_u = MatchedDisturbance(t);
    w_x = UnmatchedDisturbance(t);

    linear = Gravity(t,x);
    xdot[1:2] = linear[1:2] + w_x;
    xdot[3:4] = linear[3:4] + u + w_u;
end