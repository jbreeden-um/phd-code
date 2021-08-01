function [u, flag] = CalculateU(t, x)
persistent a_max
global outdata is_error sim_case
if isempty(a_max)
    file = load('InData/SimData.mat');
    a_max = file.a_max;
end

% CBF Conditions
if sim_case == 1
    [A, b] = Constraint_AvoidanceInfNorm(x);
elseif sim_case == 2
    [A, b] = Constraint_Avoidance(x);
elseif sim_case == 3
    [A, b] = Constraint_AvoidanceNoLim(x);
else
    error('Invalid sim_case');
end
if ~isempty(A)
    for i=1:size(A,1)
        scale = eps/max(abs(A(i,:)))*1e12;
        if ~isinf(scale)
            A(i,:) = A(i,:)*scale;
            b(i) = b(i)*scale;
        end
    end
end

% CLF Condition
r = x(1:3);
v = x(4:6);
[sdot, rg] = GroundMotion(x);
kp = 0.5;
kd = 0.1;
vd = -kp*(r - rg);
V = 1/2*dot(r-rg, r-rg) + 1/2*kd*dot(v-vd, v-vd);
Lf = (r-rg)'*v;
Lg = kd*(v-vd)';
alpha = -0.1*V;
b_clf = alpha - Lf;

J = eye(4);
J(4,4) = 10;
F = zeros(4,1);
A_cond = [A, zeros(size(A,1), 1); Lg, 1];
b_cond = [b; b_clf];
lower = -[a_max; a_max; a_max; Inf];
upper =  [a_max; a_max; a_max; Inf];
if sim_case == 3
    linopts = optimoptions('linprog', 'Algorithm', 'dual-simplex', 'Display', 'off');
    try
        [~, ~, flagc] = linprog(zeros(4, 1), A_cond, b_cond, [], [], lower, upper, linopts);
    catch
        % This indicates conflicting constraints which have no solution
        flagc = 0;
    end
    % If there is no solution, relax the control input constraints.
    % We only do this in case 3, because in all the other cases, a feasible control input
    % is guaranteed to exist.
    if flagc~= 1
        lower = lower*Inf;
        upper = upper*Inf;
    end
end


opts = optimset('Display', 'off', 'MaxIter', 2000, 'TolCon', 1e-8, 'ScaleProblem', 'obj-and-constr');
[out, ~, flag] = quadprog(J, F, A_cond, b_cond, [], [], lower, upper, [], opts);

if flag < 0
    is_error = 1;

    linopts = optimoptions('linprog', 'Algorithm', 'dual-simplex', 'Display', 'off');
    try
        [out, ~, flagc] = linprog(zeros(4, 1), A_cond, b_cond, [], [], lower, upper, linopts);
    catch
        % This indicates conflicting constraints which have no solution
        flagc = 0;
    end
    if flagc==1
        % No worries. quadprog is just having a bad day.
        disp(['Solution exists at t = ', num2str(t), ', but quadprog cannot find it. Using linprog instead. Flag was ', num2str(flag)]);
    else
        % There is no solution. Choose the best solution available.
        [u, ~, ~] = linprog(sum(A,1), [], [], [], [], -a_max*[1;1;1], a_max*[1;1;1], linopts);
        out = [u; 0];
        disp(['No valid input exists at t = ', num2str(t), ', linprog failure code was ' num2str(flag), '. Instead minimizing the constraints.']);
    end

end

u = [out(1:3); sdot];

outdata.u = [u; out(4)];

end