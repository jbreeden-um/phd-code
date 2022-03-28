function x2 = UpdateX(t1, t2, x1, u)
f = [0,1,0,0;
    0,0,0,0;
    0,0,0,1;
    0,0,0,0];
g = [0,0;1,0;0,0;0,1];
x2 = x1 + (f*x1 + g*u)*(t2-t1);
end