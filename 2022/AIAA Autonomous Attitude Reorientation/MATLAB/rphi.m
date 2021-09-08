function out = rphi(t,x,u)
% Function to call phi with only a 2DOF input (pointing vector instead of full quaternion)
out = phi(t,[1;0;0;0;x(4:end)],u,x(1:3));
end