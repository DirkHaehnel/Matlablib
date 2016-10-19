function Pfeil(an,en,base,pos,cax,al,arrowlen,rad)

flag = ishold;

if nargin<6 || isempty(al)
    al = pi/5;
end
if nargin<3 || isempty(base)
    base = 1;
end
if nargin<7 || isempty(arrowlen)
    arrowlen = 0.3;
else
    arrowlen = 0.3*arrowlen;
end
if nargin<4 || isempty(pos)
    pos = [0 0 0];
end
if nargin<8 || isempty(rad)
    rad = 0.05;
else
    rad = 0.05*rad;
end
if nargin<5 || isempty(cax)
    cax = [0 1];
end

len = sqrt(sum((en - an).^2));
dir = (en - an)/len;
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
dir = len*arrowlen*dir;
u = 0:pi/50:2*pi; 
t = [0:0.1:1]'; col = ones(size(t));

en = an + base*(en-an);

surf(pos(1) + an(1) + col*(nx(1)*cos(u) + ny(1)*sin(u)) + t*ones(size(u))*(en(1) - an(1) - dir(1)),...
    pos(2) + an(2) + col*(nx(2)*cos(u) + ny(2)*sin(u)) + t*ones(size(u))*(en(2) - an(2) - dir(2)),...
    pos(3) + an(3) + col*(nx(3)*cos(u) + ny(3)*sin(u)) + t*ones(size(u))*(en(3) - an(3) - dir(3)),...
    cax(1) + diff(cax)*(1-arrowlen)*t*ones(size(u)));

hold on

surf(pos(1) + an(1) + [1;1+tan(al)]*(nx(1)*cos(u) + ny(1)*sin(u)) + [1;1]*ones(size(u))*(en(1) - an(1) - dir(1)),...
    pos(2) + an(2) + [1;1+tan(al)]*(nx(2)*cos(u) + ny(2)*sin(u)) + [1;1]*ones(size(u))*(en(2) - an(2) - dir(2)),...
    pos(3) + an(3) + [1;1+tan(al)]*(nx(3)*cos(u) + ny(3)*sin(u)) + [1;1]*ones(size(u))*(en(3) - an(3) - dir(3)),...
    cax(1) + diff(cax)*[1-arrowlen;1-arrowlen]*ones(size(u)));

surf(pos(1) + an(1) + [0;1]*(nx(1)*cos(u) + ny(1)*sin(u)),...
    pos(2) + an(2) + [0;1]*(nx(2)*cos(u) + ny(2)*sin(u)),...
    pos(3) + an(3) + [0;1]*(nx(3)*cos(u) + ny(3)*sin(u)),...
    cax(1)*[1;1]*ones(size(u)));

surf(pos(1) + en(1) - dir(1) + (1+tan(al))*t*(nx(1)*cos(u) + ny(1)*sin(u)) + (1 - t)*ones(size(u))*dir(1),...
    pos(2) + en(2) - dir(2) + (1+tan(al))*t*(nx(2)*cos(u) + ny(2)*sin(u)) + (1 - t)*ones(size(u))*dir(2),...
    pos(3) + en(3) - dir(3) + (1+tan(al))*t*(nx(3)*cos(u) + ny(3)*sin(u)) + (1 - t)*ones(size(u))*dir(3),...
    cax(1) + diff(cax)*(1-arrowlen+arrowlen*(1-t))*ones(size(u)));

% axis image
% shading interp
% axis off
% colormap hot
% caxis(sort(cax));
% camlight

if ~flag hold off; end