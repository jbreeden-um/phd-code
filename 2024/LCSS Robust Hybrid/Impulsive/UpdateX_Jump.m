function x = UpdateX_Jump(t, x0, u, use_disturbance)
x = x0;
if nargin==4 && use_disturbance
    x(3:4) = x(3:4) + u + disturbance_model(t, u);
else
    x(3:4) = x(3:4) + u;
end
end

function out = disturbance_model(t, u)
persistent generator
if t==0
    generator = RandStream('mt19937ar','Seed',1);
end
direction = generator.randn(2,1);
direction = direction / norm(direction);
out = 0.05*generator.rand*norm(u)*direction;
end