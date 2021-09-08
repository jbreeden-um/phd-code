function v = RotQ(u, q)
%RotQ Rotates a vector u by q
%    Quaternions should be written with the scalar element first
%    Both vectors and quaternions should be written as column arrays
n = size(u);
qu = [zeros(1,n(2)); u];
if ~contains(class(q), 'sym')
    v = QxQ(QxQ(q, qu), QInv(q));
else
    v = QxQ(QxQ(q, qu), QConj(q));
end
v = v(2:4,:);
end

% function v = RotQ(u, q)
% %RotQ Rotates a vector u by q
% %    Quaternions should be written with the scalar element first
% %    Both vectors and quaternions should be written as column arrays
% qu = [0; u];
% v = QxQ(QxQ(q, qu), QInv(q));
% v(1)
% v = v(2:4);
% end