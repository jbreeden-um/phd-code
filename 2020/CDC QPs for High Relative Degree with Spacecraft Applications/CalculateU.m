function [u, flag] = CalculateU(t, x, Constraints)

n = length(Constraints);
A = zeros(n, 6);
b = zeros(n, 1);
for index=1:n
    C = Constraints(index);
    if C.t1-C.tau <= t && t < C.t2
        if t < C.t1
            type = 'future';
        else
            type = 'always';
        end
        
        if isequal(Constraints(index).mode, 'proximity')
            [Ai, bi] = Constraint_Proximity(x, C, t, type);
        elseif isequal(Constraints(index).mode, 'boresight')
            [Ai, bi] = Constraint_Boresight(x, C, t, type);
        elseif isequal(Constraints(index).mode, 'asteroid')
            [Ai, bi] = Constraint_Asteroid(x, C, t, type);
        else
            error('I do not recognize that type of constraint');
        end
    
        if bi==inf
            Ai = zeros(size(Ai));
            bi = 0;
        elseif bi==-inf
            error('Constraint is not achievable');
        end

        A(index,:) = Ai;
        b(index) = bi;
    end
end

% I put the weights on the delta variables here instead of inside the cost function to 
% keep the cost function better conditioned.
A(1:(n-2), (7):(6+n-2)) = 1e-6*eye(n-2); % Relaxation variables

omega = x(10:12);
J = diag([1000, 1000, 1000, 1, 1, 1, 1000*ones(1,n-2)]);
gamma = 2;
F = [zeros(3,1); gamma/2*omega; zeros(n-2, 1)];

opts = optimset('Display', 'off', 'MaxIter', 500);

[u, ~, flag] = quadprog(J, F, A, b, [], [], [], [], [], opts);

if flag < 1
    linopts = optimoptions('linprog', 'Algorithm', 'dual-simplex', 'Display', 'off');
    try
        [~, ~, flagc] = linprog(zeros(length(J), 1), A, b, [], [], [], [], linopts);
    catch
        % This indicates conflicting constraints which have no solution
        flagc = 0;
    end
    if flagc==1
        % No worries. quadprog is just having a bad day.
    else
        disp(['No possible input at t = ', num2str(t), ', failure code was ' num2str(flag)]);
    end
end

end