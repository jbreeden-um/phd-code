offset = 0;

time = load('time.42');
Perturb = load('EnvTrq00.42');
Safety = load('Safety.42');
AttErr = load('AttErr.42');
qbn = load('qbn.42');
wbn = load('wbn.42');
Wheels = load('Whl.42');

time = time(1:(end-offset),:);
Perturb = Perturb(1:(end-offset),:);
Safety = Safety(1:(end-offset),:);
AttErr = AttErr(1:(end-offset),:);
qbn = qbn(1:(end-offset),:);
wbn = wbn(1:(end-offset),:);
Wheels = Wheels(1:(end-offset),:);

figure(1); clf;
plot(time, qbn);
xlabel 'Time (s)'; ylabel 'qbn';

figure(2); clf;
plot(time, rad2deg(AttErr));
xlabel 'Time (s)'; ylabel 'Attitude Error B-Frame (deg)';

figure(3); clf;
plot(time, Wheels(:,[3,6,9,12])*1e3); hold on;
plot(time, ones(size(time)).*0.7, 'r--')
plot(time, -ones(size(time)).*0.7, 'r--')
xlabel 'Time (s)'; ylabel 'Commanded Wheel Torques (mNm)';

figure(4); clf;
plot(time, Perturb(:,7:9)); hold on;
plot(time, vecnorm(Perturb(:,7:9),2,2), 'LineWidth', 2);
xlabel 'Time (s)'; ylabel 'Disturbance Torque B-Frame (N)';
% We only actually baselined for a disturbance torque of up to 7.19e-6 Nm.
% The above plot had torques as large as 9.02e-6 Nm.

mu = 0.00172;
h = Safety(:,[1,3,5]);
hdot = Safety(:,[2,4,6]);
H = h + absSq(hdot)/(2*mu);
figure(5); clf;
plot(time, h(:,1)); hold on;
plot(time, H(:,1), '--');
xlabel 'Time (s)';
ylabel 'Instrument Safety'

figure(6); clf;
plot(time, h(:,2)); hold on;
plot(time, H(:,2), '--');
xlabel 'Time (s)';
ylabel 'ST Sun Safety'

figure(7); clf;
plot(time, h(:,3)); hold on;
plot(time, H(:,3), '--');
xlabel 'Time (s)';
ylabel 'ST Lunar Safety'

figure(8); clf;
plot(time, Safety(:,7));
xlabel 'Time (s)'; ylabel 'Energy Safety';

figure(9); clf;
plot(time, rad2deg(wbn)); hold on;
plot(time, rad2deg(vecnorm(wbn,2,2)));
xlabel 'Time (s)'; ylabel 'Angular Rate (deg/s)';
% FYI the rate cutoff is 1 deg/s on the largest axis and 1.65 deg/s on the smallest axis