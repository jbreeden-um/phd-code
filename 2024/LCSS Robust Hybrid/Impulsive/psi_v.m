function v = psi_v(tau, t, x, tolerance, x0)
% The convergence law implemented in the following file is identical to the one used in
% the prior 2023 paper.

if nargin==4
    % The fifth argument x0 is to ensure that the number of MPC steps used is the same in
    % every call to fmincon. If psi_v becomes nonsmooth, then fmincon gets upset.
    x0 = x;
end
persistent P
if t==0
    P = compute_V(0, [], 1);
end
global psi_v_Nmax

dt = tau - t;
if dt < 0
    error('psi_v is only valid for future times');
end
if sum(isnan(x))
    error('x is NaN in psi_v');
end

mu = 398600e9;
Vd = @(x1, xc) (xc(3) - x1(3))*(P(1,1)*(xc(1) - x1(1)) + P(1,2)*(xc(2) - x1(2)) + P(1,3)*(xc(3) - x1(3)) + P(1,4)*(xc(4) - x1(4))) - (xc(2) - x1(2))*(P(2,3)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + P(2,4)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2)) - P(1,2)*(xc(3) - x1(3)) - P(2,2)*(xc(4) - x1(4))) - (xc(3) - x1(3))*(P(3,3)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + P(3,4)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2)) - P(1,3)*(xc(3) - x1(3)) - P(2,3)*(xc(4) - x1(4))) - (xc(4) - x1(4))*(P(3,4)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + P(4,4)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2)) - P(1,4)*(xc(3) - x1(3)) - P(2,4)*(xc(4) - x1(4))) - ((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2))*(P(1,3)*(xc(1) - x1(1)) + P(2,3)*(xc(2) - x1(2)) + P(3,3)*(xc(3) - x1(3)) + P(3,4)*(xc(4) - x1(4))) - ((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2))*(P(1,4)*(xc(1) - x1(1)) + P(2,4)*(xc(2) - x1(2)) + P(3,4)*(xc(3) - x1(3)) + P(4,4)*(xc(4) - x1(4))) - (xc(1) - x1(1))*(P(1,3)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + P(1,4)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2)) - P(1,1)*(xc(3) - x1(3)) - P(1,2)*(xc(4) - x1(4))) + (xc(4) - x1(4))*(P(1,2)*(xc(1) - x1(1)) + P(2,2)*(xc(2) - x1(2)) + P(2,3)*(xc(3) - x1(3)) + P(2,4)*(xc(4) - x1(4)));
Vdd = @(x1, xc) 2*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2))*(P(3,3)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + P(3,4)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2)) - P(1,3)*(xc(3) - x1(3)) - P(2,3)*(xc(4) - x1(4))) - ((mu*xc(4))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(4))/(x1(1)^2 + x1(2)^2)^(3/2) - (3*mu*xc(2)*(2*xc(1)*xc(3) + 2*xc(2)*xc(4)))/(2*(xc(1)^2 + xc(2)^2)^(5/2)) + (3*mu*x1(2)*(2*x1(1)*x1(3) + 2*x1(2)*x1(4)))/(2*(x1(1)^2 + x1(2)^2)^(5/2)))*(P(1,4)*(xc(1) - x1(1)) + P(2,4)*(xc(2) - x1(2)) + P(3,4)*(xc(3) - x1(3)) + P(4,4)*(xc(4) - x1(4))) - ((mu*xc(3))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(3))/(x1(1)^2 + x1(2)^2)^(3/2) - (3*mu*xc(1)*(2*xc(1)*xc(3) + 2*xc(2)*xc(4)))/(2*(xc(1)^2 + xc(2)^2)^(5/2)) + (3*mu*x1(1)*(2*x1(1)*x1(3) + 2*x1(2)*x1(4)))/(2*(x1(1)^2 + x1(2)^2)^(5/2)))*(P(1,3)*(xc(1) - x1(1)) + P(2,3)*(xc(2) - x1(2)) + P(3,3)*(xc(3) - x1(3)) + P(3,4)*(xc(4) - x1(4))) + 2*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2))*(P(3,4)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + P(4,4)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2)) - P(1,4)*(xc(3) - x1(3)) - P(2,4)*(xc(4) - x1(4))) - 2*(xc(3) - x1(3))*(P(1,3)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + P(1,4)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2)) - P(1,1)*(xc(3) - x1(3)) - P(1,2)*(xc(4) - x1(4))) - 2*(xc(4) - x1(4))*(P(2,3)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + P(2,4)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2)) - P(1,2)*(xc(3) - x1(3)) - P(2,2)*(xc(4) - x1(4))) - ((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2))*(P(1,1)*(xc(1) - x1(1)) + P(1,2)*(xc(2) - x1(2)) + P(1,3)*(xc(3) - x1(3)) + P(1,4)*(xc(4) - x1(4))) - ((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2))*(P(1,2)*(xc(1) - x1(1)) + P(2,2)*(xc(2) - x1(2)) + P(2,3)*(xc(3) - x1(3)) + P(2,4)*(xc(4) - x1(4))) - (xc(1) - x1(1))*(P(1,3)*((mu*xc(3))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(3))/(x1(1)^2 + x1(2)^2)^(3/2) - (3*mu*xc(1)*(2*xc(1)*xc(3) + 2*xc(2)*xc(4)))/(2*(xc(1)^2 + xc(2)^2)^(5/2)) + (3*mu*x1(1)*(2*x1(1)*x1(3) + 2*x1(2)*x1(4)))/(2*(x1(1)^2 + x1(2)^2)^(5/2))) + P(1,4)*((mu*xc(4))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(4))/(x1(1)^2 + x1(2)^2)^(3/2) - (3*mu*xc(2)*(2*xc(1)*xc(3) + 2*xc(2)*xc(4)))/(2*(xc(1)^2 + xc(2)^2)^(5/2)) + (3*mu*x1(2)*(2*x1(1)*x1(3) + 2*x1(2)*x1(4)))/(2*(x1(1)^2 + x1(2)^2)^(5/2))) + P(1,1)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + P(1,2)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2))) - (xc(2) - x1(2))*(P(2,3)*((mu*xc(3))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(3))/(x1(1)^2 + x1(2)^2)^(3/2) - (3*mu*xc(1)*(2*xc(1)*xc(3) + 2*xc(2)*xc(4)))/(2*(xc(1)^2 + xc(2)^2)^(5/2)) + (3*mu*x1(1)*(2*x1(1)*x1(3) + 2*x1(2)*x1(4)))/(2*(x1(1)^2 + x1(2)^2)^(5/2))) + P(2,4)*((mu*xc(4))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(4))/(x1(1)^2 + x1(2)^2)^(3/2) - (3*mu*xc(2)*(2*xc(1)*xc(3) + 2*xc(2)*xc(4)))/(2*(xc(1)^2 + xc(2)^2)^(5/2)) + (3*mu*x1(2)*(2*x1(1)*x1(3) + 2*x1(2)*x1(4)))/(2*(x1(1)^2 + x1(2)^2)^(5/2))) + P(1,2)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + P(2,2)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2))) - (xc(3) - x1(3))*(P(3,3)*((mu*xc(3))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(3))/(x1(1)^2 + x1(2)^2)^(3/2) - (3*mu*xc(1)*(2*xc(1)*xc(3) + 2*xc(2)*xc(4)))/(2*(xc(1)^2 + xc(2)^2)^(5/2)) + (3*mu*x1(1)*(2*x1(1)*x1(3) + 2*x1(2)*x1(4)))/(2*(x1(1)^2 + x1(2)^2)^(5/2))) + P(3,4)*((mu*xc(4))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(4))/(x1(1)^2 + x1(2)^2)^(3/2) - (3*mu*xc(2)*(2*xc(1)*xc(3) + 2*xc(2)*xc(4)))/(2*(xc(1)^2 + xc(2)^2)^(5/2)) + (3*mu*x1(2)*(2*x1(1)*x1(3) + 2*x1(2)*x1(4)))/(2*(x1(1)^2 + x1(2)^2)^(5/2))) + P(1,3)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + P(2,3)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2))) - (xc(4) - x1(4))*(P(3,4)*((mu*xc(3))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(3))/(x1(1)^2 + x1(2)^2)^(3/2) - (3*mu*xc(1)*(2*xc(1)*xc(3) + 2*xc(2)*xc(4)))/(2*(xc(1)^2 + xc(2)^2)^(5/2)) + (3*mu*x1(1)*(2*x1(1)*x1(3) + 2*x1(2)*x1(4)))/(2*(x1(1)^2 + x1(2)^2)^(5/2))) + P(4,4)*((mu*xc(4))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(4))/(x1(1)^2 + x1(2)^2)^(3/2) - (3*mu*xc(2)*(2*xc(1)*xc(3) + 2*xc(2)*xc(4)))/(2*(xc(1)^2 + xc(2)^2)^(5/2)) + (3*mu*x1(2)*(2*x1(1)*x1(3) + 2*x1(2)*x1(4)))/(2*(x1(1)^2 + x1(2)^2)^(5/2))) + P(1,4)*((mu*xc(1))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(1))/(x1(1)^2 + x1(2)^2)^(3/2)) + P(2,4)*((mu*xc(2))/(xc(1)^2 + xc(2)^2)^(3/2) - (mu*x1(2))/(x1(1)^2 + x1(2)^2)^(3/2)));
Vddd_max = 4.11e-6*compute_V(t,x0);
margin = Vddd_max*dt^2/2;

if margin < tolerance
    % If the margin is sufficiently small, use the explicit formula
    v = Vd(x, get_center(t)) + max(0, Vdd(x, get_center(t)))*dt + margin;
else
    % If the margin would be too big, use one-step MPC
    ts = sqrt(tolerance/Vddd_max*2); % the discretization margin required to achieve the desired level of margin
    Ns = ceil(dt/ts);
    
    % Because safety is prioritized over convergence and we require the controller to run
    % sufficiently fast, we enforce a maximum number of MPC predictions for the
    % stabilizing part of the control law.
    % disp(['Desired steps = ' num2str(Ns), ', Allowed steps = ' num2str(psi_v_Nmax)]);
    Ns = min(Ns, psi_v_Nmax);
    
    tspan = linspace(t, tau, Ns+1);
    
    % The actual margin that would result from the number of samples we applied
    delta = tspan(2) - tspan(1);
    margin = Vddd_max*delta^2/2;
    
    % Relaxation of convergence tolerance (part of the term d in the controller)
    margin = min(margin, tolerance);
    
    [~, x_sim] = ode45(@two_body, tspan, x);
    Vdot = zeros(Ns,1);
    for i=1:Ns
        Vdot(i) = Vd(x_sim(i,:),get_center(tspan(i))) + max(0, Vdd(x_sim(i,:),get_center(tspan(i))))*delta;
    end
    if nargin==4
        v = max(Vdot) + margin;
    else
        v = Vdot + margin;
    end
end
end