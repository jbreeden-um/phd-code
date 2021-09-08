function x = linprog_exact(F, A, b, lower, upper)
F = F(:)';
n = length(F);
m = size(A,1) + size(lower,1);
if m < n
    error('This function requires more constraints than variables');
end
if size(lower,1) ~= size(upper,1)
    error('This function requires the same number of upper and lower bounds')
end
if sum(isinf(lower))+sum(isinf(upper)) > 0
    error('This function does not allow Inf bounds. Please use A*x <= b instead');
end
if size(lower,1) ~= n
    warning('This function is only guaranteed accurate when upper and lower bounds are used');
end

Al = -eye(n);
bl = -lower;
Au = eye(n);
bu = upper;
A_all = [A; Al; Au];
b_all = [b; bl; bu];

Ac = zeros(size(lower,1),size(lower,1));
bc = zeros(size(lower));
skip = zeros(size(lower));
for i=1:size(lower,1)
    if F(i) < 0
        Ac(i,i) = 1;
        bc(i) = upper(i);
    elseif F(i) > 0
        Ac(i,i) = -1;
        bc(i) = -lower(i);
    else
        skip(i) = 1;
    end
end
Ac = Ac(~skip,:);
bc = bc(~skip,:);

At = [A; Ac];
bt = [b; bc];

m = length(bt);
if m < n
    error('This function requires at least as many linearly independent constraints as there are variables');
end
indices = nchoosek(1:m, n);
res = Inf;
x = NaN*ones(n,1);
for i=size(indices,1)
    A_guess = At(indices(i,:),:);
    b_guess = bt(indices(i,:));
    if abs(det(A_guess)) < 1e-10
        continue;
    end
    x_guess = A_guess\b_guess;
    if sum(A_all*x_guess - b_all > 0) == 0
        res_guess = F*x_guess;
        if res_guess < res
            res = res_guess;
            x = x_guess;
        end
    end
end

end