
input = 'GradLgH.txt';


fid = fopen(input);
line = fgetl(fid);
expression = '';
while ischar(line)
    expression = [expression, strip(line, 'left')];
    line = fgetl(fid);
end

newexpression = replace(expression, 'rx[t]', 'r(1)');
newexpression = replace(newexpression, 'ry[t]', 'r(2)');
newexpression = replace(newexpression, 'rz[t]', 'r(3)');
newexpression = replace(newexpression, 'sx', 's(1)');
newexpression = replace(newexpression, 'sy', 's(2)');
newexpression = replace(newexpression, 'sz', 's(3)');
newexpression = replace(newexpression, 'ux', 'u(1)');
newexpression = replace(newexpression, 'uy', 'u(2)');
newexpression = replace(newexpression, 'uz', 'u(3)');
newexpression = replace(newexpression, 'wx[t]', 'w(1)');
newexpression = replace(newexpression, 'wy[t]', 'w(2)');
newexpression = replace(newexpression, 'wz[t]', 'w(3)');
% newexpression = replace(newexpression, 'ArcTan[r(1), r(2)]', 'atan2(r(2), r(1))');
% newexpression = replace(newexpression, 'ArcTan[-r(1), -r(2)]', 'atan2(-r(2), -r(1))');
% newexpression = regexprep(newexpression, 'Cos\[\w+\]', 'cos(${$0(5:end-1)})');
% newexpression = regexprep(newexpression, 'Sin\[\w+\]', 'sin(${$0(5:end-1)})');
newexpression = replace(newexpression, ' + ', '+');
newexpression = replace(newexpression, ' - ', '-');

newexpression = replace(newexpression, 'sign[', 'sign(');
newexpression = replace(newexpression, 'abs[', 'abs(');
newexpression = replace(newexpression, ']', ')');
newexpression = replace(newexpression, '{', '[');
newexpression = replace(newexpression, '}', '];');
newexpression = replace(newexpression, ', ', ',');
newexpression = replace(newexpression, ' ', '*');
newexpression = replace(newexpression, ')(', ')*(');
newexpression = replace(newexpression, ')sign', ')*sign');
newexpression = replace(newexpression, ')abs', ')*abs');


disp(newexpression)

fclose(fid);