function [a, b, c] = EndSimEarly(t, y)
persistent var
if isempty(var)
    var = 1;
end
a = var;
b = 1;
c = 0;
end

% To end the simulation early, set a breakpoint here, and set var equal to 0 or -1.