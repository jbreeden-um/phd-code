
include("Gravity.jl");

function UpdateX!(xdot, x, p, t)
    Controller = p[1];
    Gravity = p[2];
    MatchedDisturbance = p[3];
    UnmatchedDisturbance = p[4];
    
    u = Controller(t,x);
    w_u = MatchedDisturbance(t);
    w_x = UnmatchedDisturbance(t);
    # @show w_u[3], w_x[3], u[3], Gravity(t,x[1:3])[3]

    xdot[1:3] = x[4:6] + w_x[1:3];
    xdot[4:6] = Gravity(t, x[1:3]) + w_u + w_x[4:6] + u;
end
