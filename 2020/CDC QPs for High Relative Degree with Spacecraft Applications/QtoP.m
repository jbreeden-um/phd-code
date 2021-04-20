function p = QtoP(q)
n = size(q)*[0; 1];
p = zeros(3,n);
for i=1:n
    a = q(:,i);
    a = a/norm(a);
    p(:,i) = a(2:4)/(1+a(1));
end
end