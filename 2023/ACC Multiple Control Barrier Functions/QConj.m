function q2 = QConj(q1)
%QConj Returns the conjugate of q1
%    Quaternions should be written with the scalar element last
q2 = [-1; -1; -1; 1].*q1;
end