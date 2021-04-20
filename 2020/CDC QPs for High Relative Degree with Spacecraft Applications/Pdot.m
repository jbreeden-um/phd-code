function d = Pdot(x, omega)
l = size(x)*[1; 0];
n = size(x)*[0; 1];
d = zeros(3, n);
if l==3
    for i=1:n
        p = x(:,i);
        p1 = p(1); p2 = p(2); p3 = p(3);
        w1 = omega(1,i); w2 = omega(2,i); w3 = omega(3,i);
        d(:,i) = [w1/4 - (p2*w3)/2 + (p3*w2)/2 + (p1^2*w1)/4 - (p2^2*w1)/4 - (p3^2*w1)/4 + (p1*p2*w2)/2 + (p1*p3*w3)/2;
            w2/4 + (p1*w3)/2 - (p3*w1)/2 - (p1^2*w2)/4 + (p2^2*w2)/4 - (p3^2*w2)/4 + (p1*p2*w1)/2 + (p2*p3*w3)/2;
            w3/4 - (p1*w2)/2 + (p2*w1)/2 - (p1^2*w3)/4 - (p2^2*w3)/4 + (p3^2*w3)/4 + (p1*p3*w1)/2 + (p2*p3*w2)/2];
    end
elseif l==4
    for i=1:n
        q = x(:,i);
        q = q/norm(q);
        q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
        w1 = omega(1,i); w2 = omega(2,i); w3 = omega(3,i);
        d(:,i) = [(q0*w1 - q2*w3 + q3*w2 + q0^2*w1 + q1^2*w1 - q0*q2*w3 + q0*q3*w2 + q1*q2*w2 + q1*q3*w3)/(2*(q0 + 1)^2);
            (q0*w2 + q1*w3 - q3*w1 + q0^2*w2 + q2^2*w2 + q0*q1*w3 - q0*q3*w1 + q1*q2*w1 + q2*q3*w3)/(2*(q0 + 1)^2);
            (q0*w3 - q1*w2 + q2*w1 + q0^2*w3 + q3^2*w3 - q0*q1*w2 + q0*q2*w1 + q1*q3*w1 + q2*q3*w2)/(2*(q0 + 1)^2)];
    end
else
    error('Invalid attitude representation');
end