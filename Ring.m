function Ring(a,b,rad,pos)

if nargin<3 | isempty(rad)
    rad = 0.1;
end
if nargin<4 | isempty(pos)
    pos = [0 0 0];
end
if sum(abs(a))==0 | sum(abs(b))==0 | sum(abs(cross(a,b)))==0
    a = a+b;
    len = sqrt(sum(a.^2));
    dir = a/len;
    k = sqrt(dir(1)^2 + dir(2)^2);
    if k == 0
        nx = [1, 0, 0];
        ny = [0, 1, 0];
    else 
        nx = [-(dir(3)*dir(1))/k, -(dir(3)*dir(2))/k, k];
        ny = [nx(2)*dir(3) - nx(3)*dir(2), nx(3)*dir(1) - nx(1)*dir(3), nx(1)*dir(2) - nx(2)*dir(1)];
    end
    nx = rad*nx*len;
    ny = rad*ny*len;
    u = 0:pi/50:2*pi; 
    t = [0:0.2:1]';
    col = ones(size(t));
    surf(pos(1) +  col*(nx(1)*cos(u) + ny(1)*sin(u)) + t*ones(size(u))*a(1),...
        pos(2) + col*(nx(2)*cos(u) + ny(2)*sin(u)) + t*ones(size(u))*a(2),...
        pos(3) + col*(nx(3)*cos(u) + ny(3)*sin(u)) + t*ones(size(u))*a(3),...
        t*ones(size(u)));
    hold on
    surf(pos(1) +  col*(nx(1)*cos(u) + ny(1)*sin(u)) - t*ones(size(u))*a(1),...
        pos(2) + col*(nx(2)*cos(u) + ny(2)*sin(u)) - t*ones(size(u))*a(2),...
        pos(3) + col*(nx(3)*cos(u) + ny(3)*sin(u)) - t*ones(size(u))*a(3),...
        t*ones(size(u)));
    surf(pos(1) + a(1) + [0;1]*(nx(1)*cos(u) + ny(1)*sin(u)),...
        pos(2) + a(2) + [0;1]*(nx(2)*cos(u) + ny(2)*sin(u)),...
        pos(3) + a(3) + [0;1]*(nx(3)*cos(u) + ny(3)*sin(u)),...
        [1;1]*ones(size(u)));
    surf(pos(1) - a(1) + [0;1]*(nx(1)*cos(u) + ny(1)*sin(u)),...
        pos(2) - a(2) + [0;1]*(nx(2)*cos(u) + ny(2)*sin(u)),...
        pos(3) - a(3) + [0;1]*(nx(3)*cos(u) + ny(3)*sin(u)),...
        [1;1]*ones(size(u)));
    hold off
else
    dir = cross(a,b); dir = dir/sqrt(sum(dir.^2));
    kappa = sum(a.^2-b.^2)/sum(a.*b)/2;
    if abs(kappa) > 0
        zeta1 = atan(-kappa+sqrt(kappa^2+1));
        zeta2 = atan(-kappa-sqrt(kappa^2+1));
        if sum(a.^2)*cos(zeta1)^2+sum(b.^2)*sin(zeta1)^2+sum(a.*b)*sin(2*zeta1)>sum(a.^2)*cos(zeta2)^2+sum(b.^2)*sin(zeta2)^2+sum(a.*b)*sin(2*zeta2)
            nx = a*cos(zeta1) + b*sin(zeta1);
            ny = a*cos(zeta2) + b*sin(zeta2);
        else
            ny = a*cos(zeta1) + b*sin(zeta1);
            nx = a*cos(zeta2) + b*sin(zeta2);
        end
    else
        nx = a;
        ny = b;
    end
    fac = max([sqrt(sum(nx.^2)) sqrt(sum(ny.^2))]);
    dir = fac*dir; 
    
    u = 0:pi/50:2*pi; 
    row = ones(size(u));
    v = (0:pi/100:2*pi)';
    col = ones(size(v));
    nn = cos(v)*nx+sin(v)*ny;
    nd = -sin(v)*nx+cos(v)*ny;
    nd = nd./(sqrt(sum(nd.^2,2))*[1 1 1]);
    nd = nn - (sum(nn.*nd,2)*[1 1 1]).*nd;
    nd = nd./(sqrt(sum(nd.^2,2))*[1 1 1])*fac;
    surf(pos(1) + nn(:,1)*row+nd(:,1)*rad*cos(u)+dir(1)*col*rad*sin(u),...
        pos(2) + nn(:,2)*row+nd(:,2)*rad*cos(u)+dir(2)*col*rad*sin(u),...
        pos(3) + nn(:,3)*row+nd(:,3)*rad*cos(u)+dir(3)*col*rad*sin(u),...
        v/2/pi*row);
end
axis image
shading interp
axis off
colormap jet
caxis([0 1]);
camlight
