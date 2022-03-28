function out = mu_func(t,x)
k = 1;
vdes = 12;
out = k*[vdes - x(2); vdes - x(4)];
end