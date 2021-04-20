function out = Quadratic_Always(t, t2, h, hdot)
% Returns maximum \ddot{h} such that a constraint is satisfied over the interval [t, t1]
% if h(t)=h and \dot{h}(t)=hdot
% Requires that t < t2, and that the constraint is currently satisfied
% QP should choose u to ensure that \ddot{h} \leq out

if t > t2
    error('Current time is outside the time interval');
end
if h > 0
    error('Constraint is not yet met');
end

M = 1e5;
if isinf(t2)
    if hdot > 0
        out = hdot^2/(2*h);
    else
        out = M;
    end
else
    if hdot > 0 && -2*h/hdot < t2
        out = hdot^2/(2*h);
    else
        out = (-2*h - 2*hdot*(t2 - t))/(t2-t)^2;
    end
end