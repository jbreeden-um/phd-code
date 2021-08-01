function xdot = UpdateX(t, x, u)
persistent SimDur h
if isempty(SimDur)
    file = load('InData/SimData.mat');
    SimDur = file.SimDur;
    h = waitbar(0, 'Simulating Scenario');
end
rvec = x(1:3);
vvec = x(4:6);

% Gravity Model is in meters, while simulation is in km
[a1, a2, a3] = gravitysphericalharmonic(rvec'*1e3, 'Custom', 16, {'InData/Eros_MatlabModel.mat' @load}, 'None');
accel = [a1; a2; a3]/1e3;

vvecdot = accel + u(1:3);
sdot = u(4);

xdot = [vvec; vvecdot; sdot];
waitbar(t/SimDur, h);
end