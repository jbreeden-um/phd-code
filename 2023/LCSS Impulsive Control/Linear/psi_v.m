function v = psi_v(tau, t, x, tolerance, x0)
if nargin==4
    % The fifth argument x0 is to ensure that the number of MPC steps used is the same in
    % every call to fmincon. If psi_v becomes nonsmooth, then fmincon gets upset.
    x0 = x;
end
global A_sys P psi_v_Nmax
persistent PA PAA e_max
if isempty(e_max)
    A = A_sys;
    PA = P*A + A'*P;
    PAA = PA*A + A'*PA;
    PAAA = PAA*A + A'*PAA;
    e_max = max(eig(PAAA));
end
dt = tau - t;
if dt < 0
    error('psi_v is only valid for future times');
end
if sum(isnan(x))
    error('x is NaN in psi_v');
end
mag = local_bound(x0);
v3 = e_max*mag^2;
margin = v3*dt^2/2;
if margin < tolerance
    % If the margin is sufficiently small, use the explicit formula
    v = x'*PA*x + max(0, x'*PAA*x)*dt + margin;
else
    % If the margin would be too big, use one-step MPC
    ts = sqrt(tolerance/v3*2); % the discretization margin required to achieve the desired level of margin
    Ns = ceil(dt/ts);
    
    % Because safety is prioritized over convergence and we require the controller to run
    % sufficiently fast, we enforce a maximum number of MPC predictions for the
    % stabilizing part of the control law.
    Ns = min(Ns, psi_v_Nmax);
    
    tspan = linspace(t, tau, Ns+1);
    
    % The actual margin that would result from the number of samples we applied
    delta = tspan(2) - tspan(1);
    margin = v3*delta^2/2;
    
    % Relaxation of convergence tolerance (part of the term d in the controller)
    margin = min(margin, tolerance);
    
    A = A_sys;
    [~, x_sim] = ode45(@(t,x) A*x, tspan, x);
    Vdot = zeros(Ns,1);
    for i=1:Ns
        Vdot(i) = x_sim(i,:)*PA*x_sim(i,:)' + max(0, x'*PAA*x)*delta;
    end
    if nargin==4
        v = max(Vdot) + margin;
    else
        v = Vdot + margin;
    end
end
end

% Note that we can potentially reduce margin by evaluating v3 exactly along the
% trajectory rather than assuming the worst possible upper bound, but that is not
% considered in this file.