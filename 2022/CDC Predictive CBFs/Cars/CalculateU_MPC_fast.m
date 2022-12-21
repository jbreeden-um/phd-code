function [u, H] = CalculateU_MPC(t,x)
persistent nlobj last_u nloptions yref
h = h_func(t,x);
H = h;

if t==-1
    k = 1;
    vdes = 12;

    dt = .1;
    % Horizon is the same as the horizon of Hstar
    p = round(2.5/dt);

    nx = 4;
    ny = 4;
    nu = 2;

    nlobj = nlmpc(nx, ny, nu);
    nlobj.Ts = dt;
    nlobj.PredictionHorizon = p;
    nlobj.ControlHorizon = p;

    nlobj.Model.StateFcn = @CarContinuous;
    nlobj.Jacobian.StateFcn = @CarJacobian;
    nlobj.Model.IsContinuousTime = true;
    nlobj.Model.OutputFcn = @(x,u) x;

    x0 = [-37; 10; -40; 10];
    u0 = [1.8782; 1.8782];
    validateFcns(nlobj, x0, u0)

    Q = k;
    R = 1;
    nlobj.MV = struct('Min',{-10,-10},'Max',{10,10});
    nlobj.Weights.OutputVariables = [0 1 0 1]*Q;
    nlobj.Weights.ManipulatedVariables = [1 1]*R;
    nlobj.Weights.ManipulatedVariablesRate = [0 0];

    nlobj.Optimization.CustomIneqConFcn = @myIneqConFunction;

    yref = [0, vdes, 0, vdes].*ones(p,1);
    nloptions = nlmpcmoveopt;
    nloptions.MVTarget = [0 0]; 
    last_u = [0; 0];
    
    nlobj.Optimization.SolverOptions.MaxIterations = 8;
end

[u, nloptions, info] = nlmpcmove(nlobj,x,last_u,yref,[],nloptions);

last_u = u;
end
function dxdt = CarContinuous(x, u)
f = [0,1,0,0;
    0,0,0,0;
    0,0,0,1;
    0,0,0,0];
g = [0,0;1,0;0,0;0,1];
dxdt = f*x + g*u;
end

function [f,g] = CarJacobian(x, u)
f = [0,1,0,0;
    0,0,0,0;
    0,0,0,1;
    0,0,0,0];
g = [0,0;1,0;0,0;0,1];
end

function C = myIneqConFunction(X,U,e,data)
p = size(X,1);
C = zeros(p,1);
for i=1:p
    C(i) = h_func_real(X(i,:));
end
end

function out = h_func_real(x)
rho = 2;
out = rho - norm(lane1(x(1)) - lane2(x(3)));
end