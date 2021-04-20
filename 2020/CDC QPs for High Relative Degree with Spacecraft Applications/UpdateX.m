function xdot = UpdateX(t, x, u)
persistent SimDur h
if isempty(SimDur)
    file = load('InData/SimData.mat');
    SimDur = file.SimDur;
    h = waitbar(0, 'Simulating Scenario');
end
rvec = x(1:3);
vvec = x(4:6);
mrp = x(7:9);
omega = x(10:12);

[a1, a2, a3] = gravitysphericalharmonic(rvec'*1e3, 'Custom', 16, {'InData/Eros_MatlabModel.mat' @load}, 'None');
accel = [a1; a2; a3]/1000;

vvecdot = accel + u(1:3);
mrpdot = Pdot(mrp, omega);
omegadot = u(4:6);

xdot = [vvec; vvecdot; mrpdot; omegadot];
waitbar(t/SimDur, h);
end