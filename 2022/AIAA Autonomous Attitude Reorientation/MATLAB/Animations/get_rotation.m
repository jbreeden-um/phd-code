function R = get_rotation(v0, v1)
v0 = v0/norm(v0);
v1 = v1/norm(v1);
axis = cross(v0, v1);
if norm(axis) < 1e-10
    axis = [0; 0; 0;];
else
    axis = axis/norm(axis);
end
angle = real(acos(dot(v0, v1)));

dq = [cos(angle/2); axis*sin(angle/2)];
% q0 = RtoQ(R0);
% q1 = QxQ(q0, dq);
R = QtoR(dq);
end