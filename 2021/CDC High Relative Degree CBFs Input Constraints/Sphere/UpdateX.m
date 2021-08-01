function xdot = UpdateX(t, x, u)
persistent SimDur h
if isempty(SimDur)
    file = load('InData/SimData.mat');
    SimDur = file.SimDur;
    h = waitbar(0, 'Simulating Scenario');
end
vvec = x(4:6);
vvecdot = u(1:3);
sdot = u(4);

xdot = [vvec; vvecdot; sdot];
waitbar(t/SimDur, h);
end