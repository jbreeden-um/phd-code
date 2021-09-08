function out = rhdot(t,x)
% Function to call h with only a 2DOF input (pointing vector instead of full quaternion)
out = hdot(t,[1;0;0;0;x(4:end)],x(1:3));
end