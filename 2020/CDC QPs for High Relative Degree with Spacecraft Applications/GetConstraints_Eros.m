function out = GetConstraints_Eros
rho = 10; % km
rhobar = 1; % km
beta = deg2rad(10); % radians

n = 6; % Do not change this! It will mess with the randomization seed.
tol = 2;
conA.mode = 'asteroid';
conA.t1 = -1;
conA.t2 = inf;
conA.scalar = rhobar; % 1 km above the surface
conA.vector = 'Eros';
conA.tau = 0;
conA.epsilon = tol^2 + 2*rhobar*tol;

conP.mode = 'proximity';
conP.t1 = -1;
conP.t2 = inf;
conP.scalar = 50;
conP.vector = [0; 0; 0];
conP.tau = 0;
conP.epsilon = 2*rho*tol - tol^2;

out(2*n+3) = conA;
out(2*n+4) = conP;

file = load('InData/Eros_Shape.mat');
Eros = (file.vertices)';
indices = randi([1, length(Eros)], n, 1);
times = randi([100, 10000], n, 1);
times = sort(times);
times(n+1) = 10000;

for i=1:4
    [out(2*i-1), out(2*i)] = GetTarget(Eros(:,indices(i)), times(i), (times(i)+times(i+1))/2, rho, beta);
end

[out(2*i+1), out(2*i+2)] = GetTarget(Eros(:,indices(i)+80), times(i), (times(i)+times(i+1))/2, rho, beta);

out = out([1:10, end-1:end]);
end

function [con1, con2] = GetTarget(z, t1, t2, rho, beta)
tol = rho*0.2;
con1.mode = 'proximity';
con1.t1 = t1;
con1.t2 = t2;
con1.scalar = rho;
con1.vector = z;
con1.tau = 60;
con1.epsilon = 2*rho*tol - tol^2;

tol = beta*0.2;
con2.mode = 'boresight';
con2.t1 = t1;
con2.t2 = t2;
con2.scalar = beta;
con2.vector = z;
con2.tau = 60;
con2.epsilon = cos(beta-tol) - cos(beta);
end