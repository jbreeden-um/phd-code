function varargout = PlotCone(v, theta, update)
rho = 1;
Alpha = 0.75;
p = [0;0;0];
ctheta = cos(theta);
r1 = rho*sqrt(1-ctheta^2);
width = rho*ctheta;
rset1 = (0:1:20)*r1/20;
rset2 = sqrt((10:-1:0))*r1/sqrt(10);
n = 20;
if nargout == 3
    n = 100;
end
[x1, y1, z1] = cylinder(rset1, n);
[x2, y2, z2] = cylinder(rset2, n);
z1 = z1*width;
z2 = z2*(rho-width) + width;
x = [x1; x2];
y = [y1; y2];
z = [z1; z2];

w = cross([0;0;1], v);
if norm(w)> 1e-6
    w = w./norm(w);
else
    w = [0; 1; 0];
end
b = acos(dot([0;0;1], v));
R = expm(skew(w)*b);
s = size(x);
origin = [x(:)'; y(:)'; z(:)'];
rotated = R*origin;
translated = rotated + p;
xp = reshape(translated(1,:), s);
yp = reshape(translated(2,:), s);
zp = reshape(translated(3,:), s);
if nargin==2
    if nargout==1
        h = surf(xp,yp,zp,'FaceAlpha',Alpha,'EdgeAlpha',0,'FaceColor',[1;0;0]);
        varargout{1} = h;
        % axis equal;
    elseif nargout==3
        varargout{1} = xp;
        varargout{2} = yp;
        varargout{3} = zp;
    else
        error('Invalid number of output arguments')
    end
elseif nargin==3
    set(update, 'XData', xp, 'YData', yp, 'ZData', zp);
else
    error('Invalid number of input arguments');
end
end