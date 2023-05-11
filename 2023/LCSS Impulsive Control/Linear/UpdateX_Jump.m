function x = UpdateX_Jump(x0, u)
x = x0;
x(3:4) = x(3:4) + u;
end