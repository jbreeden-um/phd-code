function [hhat,A,b,h] = CBF_icosa(t,x,rhohat,last_u)
persistent vectors
if t==0
    vectors = make_die;
end
alpha = 0.0004;
gamma = 120;
rho = 1e4;

r = x(1:3);
v = x(4:6);

x0 = get_center(t);
r0 = x0(1:3);
v0 = x0(4:6);
u0 = [zeros(3), eye(3)]*two_body(t,x0);
u_unc = [zeros(3), eye(3)]*two_body(t,x);
kappa = vectors'*(r-r0) - rho;
kappa_dot = vectors'*(v-v0);
kappa_ddot_unc = vectors'*(u_unc - u0);

h = kappa + gamma*kappa_dot;

l_hr = 1;
l_hv = gamma;

global wg_rate wc_max
wg = wg_rate*norm(last_u); % this is assuming continuous control so as to avoid writing a QCQP
w = wg + wc_max;
mu = 398600e9;
a = 42164e3;
l_fr = 2*mu/a^3;
l_fv = 0;

hhat = h + l_hr*rhohat(1) + l_hv*rhohat(2);
A = gamma*vectors';
b = -alpha*hhat - kappa_dot - gamma*kappa_ddot_unc - l_hr*rhohat(2) ...
    - l_hv*(l_fr*rhohat(1) + l_fv*rhohat(2) + w);
end