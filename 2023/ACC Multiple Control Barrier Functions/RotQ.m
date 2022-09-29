function v = RotQ(u, q)
%RotQ Rotates a vector u by q
%    Quaternions should be written with the scalar element last
%    Both vectors and quaternions should be written as column arrays
n = size(u);
qu = [u; zeros(1,n(2))];
if ~contains(class(q), 'sym')
    v = QxQ(QxQ(q, qu), QInv(q));
else
    v = QxQ(QxQ(q, qu), QConj(q));
end
v = v(1:3,:);
end