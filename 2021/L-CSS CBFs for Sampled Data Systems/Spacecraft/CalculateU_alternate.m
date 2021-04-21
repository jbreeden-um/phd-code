% With CalculateU as in this file, the results for all the simulations that work will be the same.
% It is possible that some of the simulations with errors will break sooner. The only difference
% between this and the main CalculateU.m file is that Aw = [eye(3); -eye(3)] here.
function u = CalculateU(t, x, dt)
global outdata is_error s_target constants

% CBF Conditions
[A, b] = Constraint(x, dt, 1);
scale = max(norm(b)/norm(A), 1e-4);
A = A/scale;
b = b/scale;

% Angular velocity conditions
Aw = [eye(3); -eye(3)];
bw = [(constants.w_max - x(4:6))/dt; (constants.w_max + x(4:6))/dt];

J = eye(3);
lower = -constants.u_max;
upper =  constants.u_max;

orth = cross(s_target, x(1:3));
orth = orth/norm(orth)*acos(dot(s_target, x(1:3)));
omega = x(4:6);

kp = 0.1;
kd = 0.5;
u0 = -kp*orth - kd*omega;
F = -2*u0;

At = [A; Aw];
bt = [b; bw];

% We will solve a feasibility problem first. This considers feasibility just with respect
% to the CBF, not the omega limitation.
linopts = optimoptions('linprog', 'Algorithm', 'dual-simplex', 'Display', 'off');
try
    % If there are NaNs or a couple other issues, this line of code could fail outright,
    % but that is not very common, so this try-catch is probably overkill.
    [lin_out, ~, lin_flag] = linprog(F, A, b, [], [], lower, upper, linopts);
catch
    lin_flag = -100;
end
if lin_flag ~= 1
    % If this code is executed, then the problem is not feasible. This will probably
    % happen a lot when method=='m0'. The reason for this is that we assumed feasibility
    % of the ZOH-CBF condition, which obviously is not a great assumption when the margin
    % is very large.
    if lin_flag == -2
        disp(['Larger u may be required']); 
            % This is a dummy warning. You do not need to do anything if this happens.
            % This is just to let you know that your choice of mu and umax is not very
            % good, but the simulation will still run.
        lower = -[Inf; Inf; Inf];
        upper =  [Inf; Inf; Inf];
    else
        disp(['No valid input exists at t = ', num2str(t), ...
            ', failure code was ' num2str(lin_flag), '. Returning zeros and skipping this step.']);
        u = [0; 0; 0];
        outdata.u = [nan; nan; nan];
        if ~is_error % only update the error code if there is not one already
            is_error = 1;
        end
        return;
    end
end
    
Aq = [At; eye(3)];
lq = [-inf*ones(size(bt)); lower];
uq = [bt; upper];

solver = osqp;
settings = solver.default_settings();
settings.verbose = 0;
solver.setup(J, F, Aq, lq, uq, settings);
result = solver.solve();
u = result.x(1:3);

% If the above methods are working as they should, then the next line should make very
% little difference. It is mostly about numerical tolerances.
% However, if lin_flag==-2, then this is essential.
u = limit(u, -constants.u_max, constants.u_max);

% Error Checking
if result.info.status_val ~= 1
    % This executes if OSQP thinks its solution is bad
    if result.info.status_val == -2
        code = 'max_iter';
    elseif result.info.status_val == -3
        code = 'primal_infeasible';
    else
        code = num2str(result.info.status_val);
    end
    e = Aq*result.x - uq;
    % Only bother fixing the error if the output was NaN. Otherwise, the code will still
    % run.
    if sum(isnan(u))
        try
            % linprog is helpful because it always obeys the constraints. If a solution
            % exists, linprog will find it.
            [lin_out, ~, lin_flag] = linprog(F, A_cond, b_cond, [], [], lower, upper, linopts);
            u = lin_out(1:3);
            final_flag = 1;
        catch
            u = [0;0;0];
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
        disp(['Unknown Error in OSQP. Code = ' num2str(result.info.status_val)]);
    end
end
outdata.u = u;

if sum(isnan(u))
    disp 'NaNs';
end

end

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