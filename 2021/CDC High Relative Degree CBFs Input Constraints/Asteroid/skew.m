function A = skew(v)
%skew Returns skew symmetric matrix of v
%   Returns the matrix by which one can premultiply an arbitrary vector to
%   obtain the cross product of the input vector and the arbitary vector
A = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end

