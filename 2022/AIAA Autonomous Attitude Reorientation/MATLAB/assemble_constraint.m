function [C, Ceq] = assemble_constraint(x, nonlin, eq)
if nargin==2
    C = nonlin(x);
    Ceq = [];
elseif nargin==3
   C = nonlin(x);
   Ceq = eq(x);
end
end