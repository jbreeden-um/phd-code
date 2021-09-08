function [x1, y1, z1] = rotate_coordinates(R, x0, y0, z0)
s = size(x0);
w0 = [x0(:), y0(:), z0(:)];
w1 = w0*R';
x1 = reshape(w1(:,1), s);
y1 = reshape(w1(:,2), s);
z1 = reshape(w1(:,3), s);
end