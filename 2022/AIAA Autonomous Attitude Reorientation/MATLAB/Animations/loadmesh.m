function handle = loadmesh
mtlname = 'CubeSatMTL.mtl';

mtl = {};
mtlcolors = [];
fid = fopen(mtlname);
line = fgetl(fid);
while ischar(line)
    x = strfind(line, 'newmtl');
    if x
        mtl{end+1} = strtrim(line((x+6):end));
    end
    x = strfind(line, 'Ka');
    if x
        if line(1)~='#'
            mtlcolors(length(mtl), 1:3) = str2num(strtrim(line((x+2):end)));
        end
    end
    line = fgetl(fid);
end

filename = 'CubeSat.obj';

nv = 1204;
nf = 2188;
v = zeros(nv, 3);
f = zeros(nf, 3);
c = zeros(nf, 3);

fid = fopen(filename,'r');
for i=1:6
    line = fgetl(fid);
end
jv = 1;
jf = 1;
current_color = [0, 0, 0];
while ischar(line)
    if isequal(line(1:2), 'v ')
        v(jv, :) = str2num(line(3:end));
        jv = jv+1;
    elseif length(line) >= 6 && isequal(line(1:6), 'usemtl')
        current_color = find_color(strtrim(line(7:end)), mtl, mtlcolors);
    elseif isequal(line(1:2), 'f ')
        line = line(3:end);
        f_parts = [0, 0, 0];
        for k=1:3
            x = strfind(line, '//');
            y = strfind(line, ' ');
            f_parts(k) = str2double(line(1:(x-1)));
            line = line((y+1):end);
        end
        f(jf, :) = f_parts;
        c(jf, :) = current_color;
        jf = jf+1;
    end
    line = fgetl(fid);
end

v = v*[1,0,0;0,-1,0;0,0,-1];
handle = trisurf(f,v(:,1)-5,v(:,2)+12,v(:,3)+18,'FaceVertexCData',c,'EdgeAlpha',0.2); axis equal;
end

function out = find_color(color_name, names, colors)
for i=1:length(names)
    if isequal(names{i}, color_name)
        out = colors(i,:);
        return;
    end
end
end