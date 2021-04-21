
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
newexpression = replace(newexpression, 'phi[t]', 'phi');
newexpression = replace(newexpression, 'ArcTan[r(1), r(2)]', 'atan2(r(2), r(1))');
newexpression = replace(newexpression, 'ArcTan[-r(1), -r(2)]', 'atan2(-r(2), -r(1))');
newexpression = regexprep(newexpression, 'Cos\[\w+\]', 'cos(${$0(5:end-1)})');
newexpression = regexprep(newexpression, 'Sin\[\w+\]', 'sin(${$0(5:end-1)})');
newexpression = replace(newexpression, ' + ', '+');
newexpression = replace(newexpression, ' - ', '-');
newexpression = replace(newexpression, '\[Pi]', 'pi');

% The sqrt is done last because it's hard to correctly identify the start and end if there
% is anything else in the expression (maybe should write a function for that)
newexpression = replace(newexpression, 'Sqrt[', 'sqrt(');
newexpression = replace(newexpression, ']', ')');
newexpression = replace(newexpression, '{', '[');
newexpression = replace(newexpression, '}', '];');
newexpression = replace(newexpression, ', ', ',');
newexpression = replace(newexpression, ' ', '*');


disp(newexpression)

fclose(fid);