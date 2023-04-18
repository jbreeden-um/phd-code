function q2 = QInv(q1)
%QInv Returns the inverse of q1
%    Quaternions should be written with the scalar element first
q2 = QConj(q1)./vecnorm(q1).^2;
end