function [u, H] = CalculateU(t,x,control_case)
if control_case==1
    [H, dHdt, dHdu] = Hstar_func(t,x);
elseif control_case==2
    [H, dHdt, dHdu] = Hecbf_func(t,x);
end

alpha = 0.01;
J = eye(3);
F = zeros(3,1);
A = dHdu;
b = -alpha*H - dHdt;
[u, fval, flag] = quadprog(J, F, A, b, [], [], [], [], [], optimset('Display', 'off'));
if flag < 1 || sum(isnan(u)) > 0
    disp(['t = ' num2str(t) ': Something went wrong in the controller']);
end

end