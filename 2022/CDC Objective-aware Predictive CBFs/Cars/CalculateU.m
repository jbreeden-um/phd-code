function [u, H] = CalculateU(t,x,control_case)
if control_case==1
    [H, dHdt, dHdu] = Hstar_func(t,x);
elseif control_case==2
    [H, dHdt, dHdu] = Hecbf_func(t,x);
end
mu = mu_func(t,x);

alpha = 1;
J = eye(2);
F = zeros(2,1);
A = dHdu;
b = -alpha*H - dHdt;
[du, fval, flag] = quadprog(J, F, A, b, [], [], [], [], [], optimset('Display', 'off'));
if flag < 1 || sum(isnan(du)) > 0
    disp(['t = ' num2str(t) ': Something went wrong in the controller']);
end

u = mu + du;

end