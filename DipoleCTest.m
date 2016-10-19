% Testing the calculation of the angular distribution of radiation of a
% dipole interacting with a cylindrical structure

close all
clear all

lamem = 570;
n0 = 1.33;
n = 1.33;
n1 = 1.5;
d0 = (1500-1e-2)/lamem*2*pi;
d = 1e-2/lamem*2*pi;
d1 = [];
drho = 0*pi;
nmax = 30;

theta = (pi/5e3:pi/2.5e3:pi)';
[fr,ff,fz] = DipoleC(theta,nmax,drho,n0,n,n1,d0,d,d1);

figure
phi = 0:pi/5e2:2*pi;
q = sqrt(n^2 - n^2*cos(theta).^2);
w = n*cos(theta);
adr = Cyl2ADR(theta,phi,fr(:,:,1:2),-nmax:nmax,n);
surf((sin(theta)*cos(phi)).*adr, (sin(theta)*sin(phi)).*adr, (cos(theta)*ones(1,length(phi))).*adr,adr);
hold on
for x = 0:pi/12:pi
    tmp = interp1(theta,adr,x);
    plot3((sin(x)*cos(phi)).*tmp, (sin(x)*sin(phi)).*tmp, (cos(x)*ones(1,length(phi))).*tmp,'linewidth',1,'color','k');
end
for x = 0:pi/12:2*pi
    tmp = interp1(phi,adr',x)';
    plot3((sin(theta)*cos(x)).*tmp, (sin(theta)*sin(x)).*tmp, cos(theta).*tmp,'linewidth',1,'color','k');
end
hold off
colormap autumn
shading interp
camlight
title('rho-dipole')
axis image
cameratoolbar

figure
adf = Cyl2ADR(theta,phi,ff(:,:,1:2),-nmax:nmax,n);
surf((sin(theta)*cos(phi)).*adf, (sin(theta)*sin(phi)).*adf, (cos(theta)*ones(1,length(phi))).*adf,adf);
hold on
for x = 0:pi/12:pi
    tmp = interp1(theta,adf,x);
    plot3((sin(x)*cos(phi)).*tmp, (sin(x)*sin(phi)).*tmp, (cos(x)*ones(1,length(phi))).*tmp,'linewidth',1,'color','k');
end
for x = 0:pi/12:2*pi
    tmp = interp1(phi,adf',x)';
    plot3((sin(theta)*cos(x)).*tmp, (sin(theta)*sin(x)).*tmp, cos(theta).*tmp,'linewidth',1,'color','k');
end
hold off
colormap autumn
shading interp
camlight
title('phi-dipole')
axis image
cameratoolbar

figure
adz = Cyl2ADR(theta,phi,fz(:,:,1:2),-nmax:nmax,n);
surf((sin(theta)*cos(phi)).*adz, (sin(theta)*sin(phi)).*adz, (cos(theta)*ones(1,length(phi))).*adz,adz);
hold on
for x = 0:pi/12:pi
    tmp = interp1(theta,adz,x);
    plot3((sin(x)*cos(phi)).*tmp, (sin(x)*sin(phi)).*tmp, (cos(x)*ones(1,length(phi))).*tmp,'linewidth',1,'color','k');
end
for x = 0:pi/12:2*pi
    tmp = interp1(phi,adz',x)';
    plot3((sin(theta)*cos(x)).*tmp, (sin(theta)*sin(x)).*tmp, cos(theta).*tmp,'linewidth',1,'color','k');
end
hold off
colormap autumn
shading interp
camlight
title('z-dipole')
axis image
cameratoolbar

return

for j=1:36
    view([j*10 30]);
    eval(['print -dpng -r300 -zbuffer adr' mint2str(j,2)]);
end
