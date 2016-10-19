function RingPfeil(a,b,al,arrowlen,pos,rad)

flag = ishold;

if nargin<3 | isempty(al)
    al = pi/5;
end
if nargin<4 | isempty(arrowlen)
    arrowlen = 0.08;
else
    arrowlen = 0.08*arrowlen;
end
if nargin<5 | isempty(pos)
    pos = [0 0 0];
end
if nargin<6 | isempty(rad)
    rad = 0.1;
else
    rad = rad*0.1;
end

if sum(abs(a))<0.1*sum(abs(b)) | sum(abs(b))<0.1*sum(abs(a)) | sum(abs(cross(a,b)))==0
    Pfeil([0 0 0],a+b,[],[],pos);    
    hold on
    Pfeil([0 0 0],-a-b,[],[],pos);
else
    dir = cross(a,b); dir = dir/sqrt(sum(dir.^2));
    nx = a; 
    ny = b - sum(nx.*b)*nx;
    fac = max([sqrt(sum(a.^2)) sqrt(sum(b.^2))]);
    dir = fac*dir; 
    
    u = 0:pi/50:2*pi; 
    row = ones(size(u));
    v = (0:pi/100:2*pi*(1-arrowlen))';
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
    
    hold on
    
    v = (2*pi*(1-arrowlen):pi/100:2*pi)';
    tmp = 1-arrowlen;
    col = ones(size(v));
    nn = cos(v)*nx+sin(v)*ny;
    nd = -sin(v)*nx+cos(v)*ny;
    nd = nd./(sqrt(sum(nd.^2,2))*[1 1 1]);
    nd = nn - (sum(nn.*nd,2)*[1 1 1]).*nd;
    nd = nd./(sqrt(sum(nd.^2,2))*[1 1 1])*fac;
    surf(pos(1) + nn(:,1)*row + nd(:,1).*rad.*(1+tan(al)).*(2*pi-v)/2/pi/arrowlen*cos(u)+...
        rad.*(1+tan(al))*dir(1)*(2*pi-v)/2/pi/arrowlen*sin(u),...
        pos(2) + nn(:,2)*row + nd(:,2).*rad.*(1+tan(al)).*(2*pi-v)/2/pi/arrowlen*cos(u)+...
        rad.*(1+tan(al))*dir(2)*(2*pi-v)/2/pi/arrowlen*sin(u),...
        pos(3) + nn(:,3)*row + nd(:,3).*rad.*(1+tan(al)).*(2*pi-v)/2/pi/arrowlen*cos(u)+...
        rad.*(1+tan(al))*dir(3)*(2*pi-v)/2/pi/arrowlen*sin(u),...
        v/2/pi*row);
    
    v = v(1);
    col = [1;1];
    nn = cos(v)*nx+sin(v)*ny;
    nd = -sin(v)*nx+cos(v)*ny;
    nd = nd./(sqrt(sum(nd.^2,2))*[1 1 1]);
    nd = nn - (sum(nn.*nd,2)*[1 1 1]).*nd;
    nd = nd./(sqrt(sum(nd.^2,2))*[1 1 1])*fac;
    surf(pos(1) + nn(:,1)+nd(:,1).*rad.*(1+tan(al)).*(2*pi-v)/2/pi/arrowlen*[0;1]*cos(u)+...
        rad.*(1+tan(al))*dir(1)*(2*pi-v)/2/pi/arrowlen*[0;1]*sin(u),...
        pos(2) + nn(:,2)+nd(:,2).*rad.*(1+tan(al)).*(2*pi-v)/2/pi/arrowlen*[0;1]*cos(u)+...
        rad.*(1+tan(al))*dir(2)*(2*pi-v)/2/pi/arrowlen*[0;1]*sin(u),...
        pos(3) + nn(:,3)+nd(:,3).*rad.*(1+tan(al)).*(2*pi-v)/2/pi/arrowlen*[0;1]*cos(u)+...
        rad.*(1+tan(al))*dir(3)*(2*pi-v)/2/pi/arrowlen*[0;1]*sin(u),...
        v/2/pi*ones(2,length(u)));
    
    axis image
    shading interp
    axis off
    colormap hot
    caxis([0 1]);
    camlight
    
    if ~flag hold off; end
end

