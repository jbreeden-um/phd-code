function [u, compute, hhat] = CalculateU(t,x,rhohat)
persistent last_u
if t==0
    last_u = [0;0;0];
end
start_u = tic;
[hhat, A, b] = CBF_icosa(t,x,rhohat,last_u);
[u, fval, flag] = quadprog(1e6*eye(3),zeros(3,1),A,b,[],[],[],[],[],optimset('Display', 'off', 'TolX', 1e-16));

last_u = u;
compute = toc(start_u);

if min(b) >= 0 && norm(u) > 1e-8
     disp(['QP solution should be zero but is not: ' num2str(norm(u))]);
end

if max(A*u - b) > 1e-6
    disp('QP solution violates constraints');
end

if min(b) < -1e7
    error('Something blew up');
end

end