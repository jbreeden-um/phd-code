function out = Quadratic_Future(t, t1, h, hdot)
% Returns maximum \ddot{h} such that a constraint is satisfied at t0 if h(t)=h and
% \dot{h}(t)=hdot
% Requires that t < t1
% QP should choose u to ensure that \ddot{h} \leq out
if t > t1
    error('This function is only to be used for future constraints');
end
out = (-2*h - 2*hdot*(t1 - t))/(t1-t)^2;
end