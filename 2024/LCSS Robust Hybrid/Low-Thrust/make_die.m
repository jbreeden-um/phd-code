function [vectors, vertices, indices] = make_die(show_plots, fig_num)
% Constructs an icosahedron (20-sided die) as our safe set
% Reference: https://en.wikipedia.org/wiki/Regular_icosahedron
if nargin==0
    show_plots = 0;
end
if nargin<=1
    fig_num = 1;
end
r = (1 + sqrt(5))/2;
radius = sqrt(1+r^2);
x = [r r -r -r];
y = [1 -1 -1 1];
z = [0 0 0 0];
vertices = [[x; y; z], [z; x; y], [y; z; x]] / radius;

if show_plots
    figure(fig_num); clf;
    plot3(vertices(1,:), vertices(2,:), vertices(3,:), '.', 'MarkerSize', 10)
    axis equal
end

indices = zeros(3, 20);
vectors = zeros(3, 20);
count = 0;
for i=1:12
    for j=(i+1):12
        for k=(j+1):12
            v1 = vertices(:,k) - vertices(:,i);
            v2 = vertices(:,j) - vertices(:,i);
            v3 = vertices(:,k) - vertices(:,j);
            if norm(v1) < 1.1 && norm(v2) < 1.1 && norm(v3) < 1.1
                count = count + 1;
                indices(:,count) = [i;j;k];
                vectors(:,count) = mean(vertices(:,[i,j,k]),2);
                if show_plots
                    hold on;
                    plot3(vertices(1,[i,j,k,i]), vertices(2,[i,j,k,i]), vertices(3,[i,j,k,i]), 'k');
                    plot3([0 vectors(1,count)], [0 vectors(2,count)], [0 vectors(3,count)], 'b', 'LineWidth', 2);
                end
            end
        end
    end
end
vectors = vectors ./ vecnorm(vectors);
end