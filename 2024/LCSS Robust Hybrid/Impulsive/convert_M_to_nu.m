function nu = convert_M_to_nu(M, e)
if M > pi, M = M - 2*pi; elseif M < -pi, M = M + 2*pi; end
% E = fzero(@(E) E-e*sin(E)-M, M);
E = my_fzero(M, e);
if E > pi, E = E - 2*pi; elseif E < -pi, E = E + 2*pi; end
if sign(E) ~= sign(M), E = -E; end
nu = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end

function x = my_fzero(M, e)
x = M;
count = 1;
while abs(x-e*sin(x)-M) > 1e-8
    x = M + e*sin(x);
    count = count+1;
    if count > 1000
        error('convert_M_to_nu did not converge');
    end
end
end