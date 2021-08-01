% Computes necessary gravitational information for the controller to use.

shape = load('Eros_Shape.mat');
scale = 1; % The maximum always occurs at lowest altitude, so this is all we have to consider.

vecs = shape.vertices*scale;

[a1, a2, a3] = gravitysphericalharmonic(vecs*1e3, 'Custom', 16, {'Eros_MatlabModel.mat' @load}, 'None');
accel_real = [a1, a2, a3]/1e3;

grav = load('Eros_MatlabModel.mat');
mu = grav.GM/1e9;
r3 = vecnorm(vecs,2,2).^3;
accel_calc = -mu./r3.*vecs;

accel_diff = accel_real - accel_calc;

fmax = max(vecnorm(accel_real,2,2))
fhatmax = max(vecnorm(accel_calc,2,2))
ftilde = max(vecnorm(accel_diff,2,2)) 
% this is already a little conservative because the safe set does not actually include the surface

save('Eros_GravApprox.mat', 'fmax', 'fhatmax', 'ftilde', 'mu');