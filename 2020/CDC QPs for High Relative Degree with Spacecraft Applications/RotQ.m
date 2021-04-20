function v = RotQ(u, q)
%RotQ Rotates a vector u by q
%    Quaternions should be written with the scalar element first
%    Both vectors and quaternions should be written as column arrays
n = size(u);
qu = [zeros(1,n(2)); u];
v = QxQ(QxQ(q, qu), QInv(q));
v = v(2:4,:);
end
