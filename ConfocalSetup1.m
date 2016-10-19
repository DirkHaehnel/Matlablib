% figure for astigmatism paper

flag = 0;

close all
NA = 1.2;
chimax = asin(NA/1.333);
chi = (0:1e2)/100*chimax;
phi = pi/2 + (0:100)/200*2*pi;
foc = 200;

red = 0.8;
lighttransparency = 0.5;

d = 2;
len = 5*foc*sin(chimax);
height = 11*foc;
surf([0 len], [-len len], -foc*ones(2,2), zeros(2,2),'facealpha', 0.1) 
hold on
surf([0 0; len len], -len*ones(2,2), -foc*[1 d; 1 d], zeros(2,2),'facealpha', 0.3) 
surf(zeros(2,2), [-len -len; len len], -foc*[1 d; 1 d], zeros(2,2),'facealpha', 0.2)
theta0 = 0.8*chimax;
rad = foc*tan(theta0);
theta = asin(1.333*sin(theta0)/1.51);

% Light cone
surf([0; -rad]*cos(phi),[0; -rad]*sin(phi),[0; -foc]*ones(1,length(phi)),0.5*ones(2,length(phi)), 'facealpha', lighttransparency)

al = pi/5;
arrowlen = 0.3;
pos = [0 0 0];
prad = 0.05;
pcol = -2;

an = [0,-rad,-foc];
en = [0,-rad,-foc]+[170 0 0];

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
nx = prad*nx*len;
ny = prad*ny*len;
dir = len*arrowlen*dir;
u = 0:pi/50:2*pi; 
t = [0:0.1:1]'; col = ones(size(t));
surf(pos(1) + an(1) + col*(nx(1)*cos(u) + ny(1)*sin(u)) + t*ones(size(u))*(en(1) - an(1) - dir(1)),...
    pos(2) + an(2) + col*(nx(2)*cos(u) + ny(2)*sin(u)) + t*ones(size(u))*(en(2) - an(2) - dir(2)),...
    pos(3) + an(3) + col*(nx(3)*cos(u) + ny(3)*sin(u)) + t*ones(size(u))*(en(3) - an(3) - dir(3)),...
    pcol*ones(size(t*u)));
surf(pos(1) + an(1) + [1;1+tan(al)]*(nx(1)*cos(u) + ny(1)*sin(u)) + [1;1]*ones(size(u))*(en(1) - an(1) - dir(1)),...
    pos(2) + an(2) + [1;1+tan(al)]*(nx(2)*cos(u) + ny(2)*sin(u)) + [1;1]*ones(size(u))*(en(2) - an(2) - dir(2)),...
    pos(3) + an(3) + [1;1+tan(al)]*(nx(3)*cos(u) + ny(3)*sin(u)) + [1;1]*ones(size(u))*(en(3) - an(3) - dir(3)),...
    pcol*ones(2,length(u)));
surf(pos(1) + an(1) + [0;1]*(nx(1)*cos(u) + ny(1)*sin(u)),...
    pos(2) + an(2) + [0;1]*(nx(2)*cos(u) + ny(2)*sin(u)),...
    pos(3) + an(3) + [0;1]*(nx(3)*cos(u) + ny(3)*sin(u)),...
    pcol*ones(2,length(u)));
surf(pos(1) + en(1) - dir(1) + (1+tan(al))*t*(nx(1)*cos(u) + ny(1)*sin(u)) + (1 - t)*ones(size(u))*dir(1),...
    pos(2) + en(2) - dir(2) + (1+tan(al))*t*(nx(2)*cos(u) + ny(2)*sin(u)) + (1 - t)*ones(size(u))*dir(2),...
    pos(3) + en(3) - dir(3) + (1+tan(al))*t*(nx(3)*cos(u) + ny(3)*sin(u)) + (1 - t)*ones(size(u))*dir(3),...
    pcol*ones(size(t*u)));

an = [0,-rad,-foc]; 
en = [0,-rad,-foc]+[170 0 0];

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
nx = prad*nx*len;
ny = prad*ny*len;
dir = len*arrowlen*dir;
u = 0:pi/50:2*pi; 
t = [0:0.1:1]'; col = ones(size(t));
surf(pos(1) + an(1) + col*(nx(1)*cos(u) + ny(1)*sin(u)) + t*ones(size(u))*(en(1) - an(1) - dir(1)),...
    pos(2) + an(2) + col*(nx(2)*cos(u) + ny(2)*sin(u)) + t*ones(size(u))*(en(2) - an(2) - dir(2)),...
    pos(3) + an(3) + col*(nx(3)*cos(u) + ny(3)*sin(u)) + t*ones(size(u))*(en(3) - an(3) - dir(3)),...
    pcol*ones(size(t*u)));
surf(pos(1) + an(1) + [1;1+tan(al)]*(nx(1)*cos(u) + ny(1)*sin(u)) + [1;1]*ones(size(u))*(en(1) - an(1) - dir(1)),...
    pos(2) + an(2) + [1;1+tan(al)]*(nx(2)*cos(u) + ny(2)*sin(u)) + [1;1]*ones(size(u))*(en(2) - an(2) - dir(2)),...
    pos(3) + an(3) + [1;1+tan(al)]*(nx(3)*cos(u) + ny(3)*sin(u)) + [1;1]*ones(size(u))*(en(3) - an(3) - dir(3)),...
    pcol*ones(2,length(u)));
surf(pos(1) + an(1) + [0;1]*(nx(1)*cos(u) + ny(1)*sin(u)),...
    pos(2) + an(2) + [0;1]*(nx(2)*cos(u) + ny(2)*sin(u)),...
    pos(3) + an(3) + [0;1]*(nx(3)*cos(u) + ny(3)*sin(u)),...
    pcol*ones(2,length(u)));
surf(pos(1) + en(1) - dir(1) + (1+tan(al))*t*(nx(1)*cos(u) + ny(1)*sin(u)) + (1 - t)*ones(size(u))*dir(1),...
    pos(2) + en(2) - dir(2) + (1+tan(al))*t*(nx(2)*cos(u) + ny(2)*sin(u)) + (1 - t)*ones(size(u))*dir(2),...
    pos(3) + en(3) - dir(3) + (1+tan(al))*t*(nx(3)*cos(u) + ny(3)*sin(u)) + (1 - t)*ones(size(u))*dir(3),...
    pcol*ones(size(t*u)));

an = [0,-rad,-foc];
en = [0,-rad,-foc]+170*[0 -foc rad]/sqrt(rad^2+foc^2);

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
nx = prad*nx*len;
ny = prad*ny*len;
dir = len*arrowlen*dir;
u = 0:pi/50:2*pi; 
t = [0:0.1:1]'; col = ones(size(t));
surf(pos(1) + an(1) + col*(nx(1)*cos(u) + ny(1)*sin(u)) + t*ones(size(u))*(en(1) - an(1) - dir(1)),...
    pos(2) + an(2) + col*(nx(2)*cos(u) + ny(2)*sin(u)) + t*ones(size(u))*(en(2) - an(2) - dir(2)),...
    pos(3) + an(3) + col*(nx(3)*cos(u) + ny(3)*sin(u)) + t*ones(size(u))*(en(3) - an(3) - dir(3)),...
    pcol*ones(size(t*u)));
surf(pos(1) + an(1) + [1;1+tan(al)]*(nx(1)*cos(u) + ny(1)*sin(u)) + [1;1]*ones(size(u))*(en(1) - an(1) - dir(1)),...
    pos(2) + an(2) + [1;1+tan(al)]*(nx(2)*cos(u) + ny(2)*sin(u)) + [1;1]*ones(size(u))*(en(2) - an(2) - dir(2)),...
    pos(3) + an(3) + [1;1+tan(al)]*(nx(3)*cos(u) + ny(3)*sin(u)) + [1;1]*ones(size(u))*(en(3) - an(3) - dir(3)),...
    pcol*ones(2,length(u)));
surf(pos(1) + an(1) + [0;1]*(nx(1)*cos(u) + ny(1)*sin(u)),...
    pos(2) + an(2) + [0;1]*(nx(2)*cos(u) + ny(2)*sin(u)),...
    pos(3) + an(3) + [0;1]*(nx(3)*cos(u) + ny(3)*sin(u)),...
    pcol*ones(2,length(u)));
surf(pos(1) + en(1) - dir(1) + (1+tan(al))*t*(nx(1)*cos(u) + ny(1)*sin(u)) + (1 - t)*ones(size(u))*dir(1),...
    pos(2) + en(2) - dir(2) + (1+tan(al))*t*(nx(2)*cos(u) + ny(2)*sin(u)) + (1 - t)*ones(size(u))*dir(2),...
    pos(3) + en(3) - dir(3) + (1+tan(al))*t*(nx(3)*cos(u) + ny(3)*sin(u)) + (1 - t)*ones(size(u))*dir(3),...
    pcol*ones(size(t*u)));

an = [0,-rad,-foc];
en = [0,-rad,-foc]+170*[0 rad foc]/sqrt(rad^2+foc^2);

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
nx = prad*nx*len;
ny = prad*ny*len;
dir = len*arrowlen*dir;
u = 0:pi/50:2*pi; 
t = [0:0.1:1]'; col = ones(size(t));
surf(pos(1) + an(1) + col*(nx(1)*cos(u) + ny(1)*sin(u)) + t*ones(size(u))*(en(1) - an(1) - dir(1)),...
    pos(2) + an(2) + col*(nx(2)*cos(u) + ny(2)*sin(u)) + t*ones(size(u))*(en(2) - an(2) - dir(2)),...
    pos(3) + an(3) + col*(nx(3)*cos(u) + ny(3)*sin(u)) + t*ones(size(u))*(en(3) - an(3) - dir(3)),...
    pcol*ones(size(t*u)));
surf(pos(1) + an(1) + [1;1+tan(al)]*(nx(1)*cos(u) + ny(1)*sin(u)) + [1;1]*ones(size(u))*(en(1) - an(1) - dir(1)),...
    pos(2) + an(2) + [1;1+tan(al)]*(nx(2)*cos(u) + ny(2)*sin(u)) + [1;1]*ones(size(u))*(en(2) - an(2) - dir(2)),...
    pos(3) + an(3) + [1;1+tan(al)]*(nx(3)*cos(u) + ny(3)*sin(u)) + [1;1]*ones(size(u))*(en(3) - an(3) - dir(3)),...
    pcol*ones(2,length(u)));
surf(pos(1) + an(1) + [0;1]*(nx(1)*cos(u) + ny(1)*sin(u)),...
    pos(2) + an(2) + [0;1]*(nx(2)*cos(u) + ny(2)*sin(u)),...
    pos(3) + an(3) + [0;1]*(nx(3)*cos(u) + ny(3)*sin(u)),...
    pcol*ones(2,length(u)));
surf(pos(1) + en(1) - dir(1) + (1+tan(al))*t*(nx(1)*cos(u) + ny(1)*sin(u)) + (1 - t)*ones(size(u))*dir(1),...
    pos(2) + en(2) - dir(2) + (1+tan(al))*t*(nx(2)*cos(u) + ny(2)*sin(u)) + (1 - t)*ones(size(u))*dir(2),...
    pos(3) + en(3) - dir(3) + (1+tan(al))*t*(nx(3)*cos(u) + ny(3)*sin(u)) + (1 - t)*ones(size(u))*dir(3),...
    pcol*ones(size(t*u)));

plot3([0 0],[-rad, -rad],-foc+[-0 120], 'linewidth', 1, 'color', [0 0 0]);
chi = (0:200)/200*theta0;
plot3(0*chi,-rad+90*sin(chi),-foc+90*cos(chi),'linewidth', 1, 'color', [0 0 0])


surf([-rad; -rad-tan(theta)*(d-1)*foc]*cos(phi),[-rad; -rad-tan(theta)*(d-1)*foc]*sin(phi),...
    [-foc; -d*foc]*ones(1,length(phi)),0.5*ones(2,length(phi)), 'facealpha', lighttransparency)
rad = rad+tan(theta)*(d-1)*foc;
plot3([0 0],[-rad, -rad],-d*foc+[-0 120], 'linewidth', 1, 'color', [0 0 0]);
chi = (0:200)/200*theta;
plot3(0*chi,-rad+90*sin(chi),-d*foc+90*cos(chi),'linewidth', 1, 'color', [0 0 0])

surf([-rad; -rad-tan(theta0)*foc]*cos(phi), [-rad; -rad-tan(theta0)*foc]*sin(phi),...
    [-d*foc; -(d+1)*foc]*ones(1,length(phi)), 0.5*ones(2,length(phi)), 'facealpha', lighttransparency)
rad = rad+tan(theta0)*foc;
plot3([0 0],[0, 0],-(d+1)*foc+[-100 700], 'linewidth', 1, 'color', [0 0 0]);
plot3([0 0],[-rad, -rad],-(d+1)*foc+[-0 120], 'linewidth', 1, 'color', [0 0 0]);
chi = (0:200)/200*theta0;
plot3(0*chi,-rad+90*sin(chi),-(d+1)*foc+90*cos(chi),'linewidth', 1, 'color', [0 0 0])


% RingPfeil
a = [0 -100 0];
b = [100 0 0];
al = pi/5;
arrowlen = 0.08;
pos = [0 0 -700];
prad = 0.1;

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
surf(pos(1) + nn(:,1)*row+nd(:,1)*prad*cos(u)+dir(1)*col*prad*sin(u),...
    pos(2) + nn(:,2)*row+nd(:,2)*prad*cos(u)+dir(2)*col*prad*sin(u),...
    pos(3) + nn(:,3)*row+nd(:,3)*prad*cos(u)+dir(3)*col*prad*sin(u),...
    pcol+v/2/pi*row);

hold on

v = (2*pi*(1-arrowlen):pi/100:2*pi)';
tmp = 1-arrowlen;
col = ones(size(v));
nn = cos(v)*nx+sin(v)*ny;
nd = -sin(v)*nx+cos(v)*ny;
nd = nd./(sqrt(sum(nd.^2,2))*[1 1 1]);
nd = nn - (sum(nn.*nd,2)*[1 1 1]).*nd;
nd = nd./(sqrt(sum(nd.^2,2))*[1 1 1])*fac;
surf(pos(1) + nn(:,1)*row + nd(:,1).*prad.*(1+tan(al)).*(2*pi-v)/2/pi/arrowlen*cos(u)+...
    prad.*(1+tan(al))*dir(1)*(2*pi-v)/2/pi/arrowlen*sin(u),...
    pos(2) + nn(:,2)*row + nd(:,2).*prad.*(1+tan(al)).*(2*pi-v)/2/pi/arrowlen*cos(u)+...
    prad.*(1+tan(al))*dir(2)*(2*pi-v)/2/pi/arrowlen*sin(u),...
    pos(3) + nn(:,3)*row + nd(:,3).*prad.*(1+tan(al)).*(2*pi-v)/2/pi/arrowlen*cos(u)+...
    prad.*(1+tan(al))*dir(3)*(2*pi-v)/2/pi/arrowlen*sin(u),...
    pcol+v/2/pi*row);

v = v(1);
col = [1;1];
nn = cos(v)*nx+sin(v)*ny;
nd = -sin(v)*nx+cos(v)*ny;
nd = nd./(sqrt(sum(nd.^2,2))*[1 1 1]);
nd = nn - (sum(nn.*nd,2)*[1 1 1]).*nd;
nd = nd./(sqrt(sum(nd.^2,2))*[1 1 1])*fac;
surf(pos(1) + nn(:,1)+nd(:,1).*prad.*(1+tan(al)).*(2*pi-v)/2/pi/arrowlen*[0;1]*cos(u)+...
    prad.*(1+tan(al))*dir(1)*(2*pi-v)/2/pi/arrowlen*[0;1]*sin(u),...
    pos(2) + nn(:,2)+nd(:,2).*prad.*(1+tan(al)).*(2*pi-v)/2/pi/arrowlen*[0;1]*cos(u)+...
    prad.*(1+tan(al))*dir(2)*(2*pi-v)/2/pi/arrowlen*[0;1]*sin(u),...
    pos(3) + nn(:,3)+nd(:,3).*prad.*(1+tan(al)).*(2*pi-v)/2/pi/arrowlen*[0;1]*cos(u)+...
    prad.*(1+tan(al))*dir(3)*(2*pi-v)/2/pi/arrowlen*[0;1]*sin(u),...
    pcol+v/2/pi*ones(2,length(u)));
    

% Objective
rad = 1.15*rad;
rho = (0:50)'/50*rad;
lens = 1.5e3;
surf(-rho*cos(phi), -rho*sin(phi), -(d+1.3)*foc-(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), 0*ones(length(rho),length(phi)), 'facealpha', 0.4)

hold off

colormap hot
caxis([-4 2])
shading interp
axis image
axis off
view([-58 11])
material shiny
camlight
% camlight

