function u = CalculateU(t, x, dt)
global outdata is_error r_target separation two_obstacles constants

u_min = constants.u_min;
u_max = constants.u_max;

% CBF Conditions
if two_obstacles
    [A1, b1] = Constraint(x, dt, 1);
    [A2, b2] = Constraint(x - [-1;1;0]*separation/sqrt(2), dt, 2);
    A = [A1; A2]; b = [b1; b2];
else
    [A, b] = Constraint(x, dt, 1);
end
scale = max(norm(b)/norm(A), 1e-4);
A = A/scale;
b = b/scale;

% CLF Condition
dx = x(1:2) - r_target;
J = eye(2);
lower = u_min;
upper = u_max;

kw = 0.5*u_max(2)/pi;
phi_des = atan2(-dx(2), -dx(1));
delta_phi = wrap_angle(phi_des - x(3));
w0 = kw*delta_phi;

kv = 0.05;
v0 = kv*norm(dx)*max(0.1, cos(delta_phi)^8);

u0 = [v0; w0];
F = -2*u0;

% We will solve a feasibility problem first
linopts = optimoptions('linprog', 'Algorithm', 'dual-simplex', 'Display', 'off');
try
    % If there are NaNs or a couple other issues, this line of code could fail outright,
    % but that is not very common, so this try-catch is probably overkill.
    [lin_out, ~, lin_flag] = linprog(F, A, b, [], [], lower, upper, linopts);
catch
    lin_flag = -100;
end
if lin_flag ~= 1
    if lin_flag == -2
        lower(1) = 0;
        [lin_out, ~, lin_flag] = linprog(F, A, b, [], [], lower, upper, linopts);
        if lin_flag==-2
            disp(['Larger u may be required']); 
                % This is a dummy warning. You do not need to do anything if the simulation works,
                % but if it doesn't work, this is a likely culprit.
            lower = -[Inf; Inf];
            upper =  [Inf; Inf];
        end
    else
        disp(['No valid input exists at t = ', num2str(t), ...
            ', failure code was ' num2str(lin_flag), '. Returning zeros and skipping this step.']);
        u = [0; 0];
        outdata.u = [nan; nan];
        if ~is_error % only update the error code if there is not one already
            is_error = 1;
        end
    return;
    end
end
    
Aq = [A; eye(2)];
lq = [-inf*ones(size(b)); lower];
uq = [b; upper];

solver = osqp;
settings = solver.default_settings();
settings.verbose = 0;
solver.setup(J, F, Aq, lq, uq, settings);
result = solver.solve();
u = result.x(1:2);

% If the above methods are working as they should, then the next line should make very
% little difference
u = limit(u, u_min, u_max);

if result.info.status_val ~= 1
    if result.info.status_val == -2
        code = 'max_iter';
    elseif result.info.status_val == -3
        code = 'primal_infeasible';
    else
        code = num2str(result.info.status_val);
    end
    e = Aq*result.x - uq;
    if sum(isnan(u))
        try
            if lin_flag==1 && 0
                u = lin_out(1:2);
            else
                [lin_out, ~, lin_flag] = linprog(F, A_cond, b_cond, [], [], lower, upper, linopts);
                u = lin_out(1:2);
            end
            final_flag = 1;
        catch
            u = [0;0];
            is_error = 100;
            final_flag = 0;
        end
        if final_flag
            disp(['Problem is feasible at t = ', num2str(t), ', but OSQP cannot find solution. ', ...
                'Code = ' code, '. Error = ' num2str(sum(e(~(e <= 0))))]); % this boolean expression is to highlight nans
        else
            disp(['Problem infeasible at t = ', num2str(t), '. ', ...
                'Code = ' code, '. Error = ' num2str(sum(e(~(e <= 0))))]);
        end
    else
        disp(['Error in OSQP. Code = ' num2str(code)]);
    end
end

outdata.u = u;

if sum(isnan(u))
    disp 'NaNs';
end

end

function out = limit(x, min, max)
out = x;
for i=1:2
    if out(i) < min(i)
        out(i) = min(i);
    elseif out(i) > max(i)
        out(i) = max(i);
    end
end
end