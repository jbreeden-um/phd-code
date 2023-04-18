function de = mpc_PE(u,R,AM,F,M,O,p1,p2,Tp,mu,dmu,eta0,dAM0,dp10,dp20,mpc)

% This function is omitted for copyright reasons.
% Interested readers should contact the authors of:
	% Gupta, R., Kalabic, U. V., Di Cairano, S., Bloch, A. M., and Kolmanovsky, I. V., 
	% "Constrained spacecraft attitude control on SO(3) using fast nonlinear model 
	% predictive control," 2015 American Control Conference, 2015, pp. 2980–2986. 
	% https://doi.org/10.1109/ACC.2015.7171188.
% We simply copied these author's code with the following changes
	% 1. The statement "-eye(3)" in the y expression was replaced with "-R_target"
	% 2. A second pointing constraint of the form of "constraint2" was added, resulting
	% in additional z and dz terms in the expression for dp1

end