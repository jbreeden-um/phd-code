% Version that actually follows the notation in the paper

%% Clean Up Workspace
clear CalculateU UpdateX Constraint_Proximity Constraint_Boresight Constraint_Asteroid Plot_Outputs
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% Set up Variables
mu = 0; % mu is no longer needed, since this simulation doesn't even call f_mu since we assume it's zero.
real_mu = 4.46176546159E-04;
ftilde_mu = 5e-5; % limit on uncertainty in gravity if we assumed no gravity
SimDur = 6000;

%% Choose Selection Method
x0_method = 'load';

if isequal(x0_method, 'generate')
    file = load('InData/RngSeed.mat');
    seed = file.seed;
    rng(seed);
    Constraints = GetConstraints_Eros;
    
    r0 = randn(3,1);
    r0 = r0/norm(r0)*12;
    rand_vec = rand(3,1);
    rand_vec = rand_vec - dot(rand_vec, r0)/dot(r0, r0)*r0;
    v0 = sqrt(real_mu/norm(r0))*rand_vec/norm(rand_vec);
    rand_vec = rand(4, 1);
    mrp0 = QtoP(rand_vec/norm(rand_vec));
    omega0 = [0; 0; 0];
elseif isequal(x0_method, 'load')
    load('InData/Constraints.mat');
end

%% Set Up Constraints
save('InData/SimData.mat', 'mu', 'Constraints', 'SimDur', 'ftilde_mu');

%% Set Up Initial Conditions
x0 = [r0; v0; mrp0; omega0];

%% Run Simulation
opts = odeset('OutputFcn', @PlotOutputs, 'Events', @EndSimEarly, 'RelTol', 1e-6, 'MaxStep', 10);
func = @(t, x) UpdateX(t, x, CalculateU(t, x, Constraints));
[t, x] = ode45(func, [0 SimDur], x0, opts);

%% Clean up
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);