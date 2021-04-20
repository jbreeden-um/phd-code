function q = PtoQ(p)
n = size(p)*[0; 1];
q = zeros(4, n);
for i=1:n
    a = p(:,i);
    s = norm(a)^2;
    rho = (1-s)/(1+s);
    v = 2*a/(1+s);
    q(:,i) = [rho; v];
end
end