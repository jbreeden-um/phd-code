function q3 = QxQ(q1, q2)
%QxQ Returns quaternion product of q1 and q2
%    Quaternions should be written with the scalar element last
if min(size(q1)) > 1
    b1 = q1(1,:);
    c1 = q1(2,:);
    d1 = q1(3,:);
    a1 = q1(4,:);
else
    b1 = q1(1);
    c1 = q1(2);
    d1 = q1(3);
    a1 = q1(4);
end
if min(size(q2)) > 1
    b2 = q2(1,:);
    c2 = q2(2,:);
    d2 = q2(3,:);
    a2 = q2(4,:);
else
    b2 = q2(1);
    c2 = q2(2);
    d2 = q2(3);
    a2 = q2(4);
end
q3 = [a1.*b2 + b1.*a2 + c1.*d2 - d1.*c2;
    a1.*c2 - b1.*d2 + c1.*a2 + d1.*b2;
    a1.*d2 + b1.*c2 - c1.*b2 + d1.*a2;
    a1.*a2 - b1.*b2 - c1.*c2 - d1.*d2];
% https://en.wikipedia.org/wiki/Quaternion
end