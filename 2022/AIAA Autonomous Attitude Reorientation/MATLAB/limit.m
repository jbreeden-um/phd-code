function out = limit(x, min, max)
out = x;
for i=1:length(x)
    if out(i) < min(i)
        out(i) = min(i);
    elseif out(i) > max(i)
        out(i) = max(i);
    end
end
end