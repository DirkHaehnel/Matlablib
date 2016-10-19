% Angular distribution of radiation of a
% dipole interacting with a cylindrical metal wire

close all
clear all

load metals
lamem = 640;
nh2o = 1.33;
nau = gold(wavelength==lamem);
metal_radius = 15;
dipole_pos_vec = 4.5;
theta = (pi/5e2:pi/2.5e2:pi)';

n0 = nau;
n = nh2o;
n1 = nh2o;
nmax = 5;

for jr = 1:length(dipole_pos_vec)
    dipole_pos = dipole_pos_vec(jr);
    d0 = metal_radius/lamem*2*pi;
    d = dipole_pos/lamem*2*pi;
    drho = d;
    d1 = [];

    [fr,ff,fz] = DipoleC(theta,nmax,drho,n0,n,n1,d0,d,d1);

    %figure
    phi = 0:pi/2.5e2:2*pi;
    q = sqrt(n^2 - n^2*cos(theta).^2);
    w = n*cos(theta);
    adr = Cyl2ADR(theta,phi,fr(:,:,1:2),-nmax:nmax,n);
    surf((sin(theta)*cos(phi)).*adr, (sin(theta)*sin(phi)).*adr, (cos(theta)*ones(1,length(phi))).*adr,adr);
    xlabel('\rho-direction'); ylabel('\phi-direction'); zlabel('\itz\rm-direction')
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
    title('\rho-dipole')
    axis image
    cameratoolbar
    ax = axis; cax = caxis; hold on; Pfeil([0 1.1*ax(3) 0],[0.5*ax(2) 1.1*ax(3) 0],[],[],cax); shading interp; hold off
    eval(['print -dpng -r300 -zbuffer adr15nmGold_' mint2str(jr) ])

    %figure
    adf = Cyl2ADR(theta,phi,ff(:,:,1:2),-nmax:nmax,n);
    surf((sin(theta)*cos(phi)).*adf, (sin(theta)*sin(phi)).*adf, (cos(theta)*ones(1,length(phi))).*adf,adf);
    xlabel('\rho-direction'); ylabel('\phi-direction'); zlabel('\itz\rm-direction')
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
    title('\phi-dipole')
    axis image
    cameratoolbar
    ax = axis; cax = caxis; hold on; Pfeil([0 1.1*ax(3) 0],[0.5*ax(2) 1.1*ax(3) 0],[],[],cax); shading interp; hold off    
    eval(['print -dpng -r300 -zbuffer adf15nmGold_' mint2str(jr) ])

    %figure
    adz = Cyl2ADR(theta,phi,fz(:,:,1:2),-nmax:nmax,n);
    surf((sin(theta)*cos(phi)).*adz, (sin(theta)*sin(phi)).*adz, (cos(theta)*ones(1,length(phi))).*adz,adz);
    xlabel('\rho-direction'); ylabel('\phi-direction'); zlabel('\itz\rm-direction')
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
    ax = axis; cax = caxis; hold on; Pfeil([0 1.1*ax(3) 0],[0.5*ax(2) 1.1*ax(3) 0],[],[],cax); shading interp; hold off
    eval(['print -dpng -r300 -zbuffer adz15nmGold_' mint2str(jr) ])

    eval(['save DipoleCyl15nmGold_' mint2str(jr) ' n0 n n1 dipole_pos_vec lamem metal_radius nmax theta phi adr adf adz fr ff fz'])
end

return

for j=1:36
    view([-j*10 30]);
    eval(['print -dpng -r300 -zbuffer adf15nmGold' mint2str(j,2)]);
end
