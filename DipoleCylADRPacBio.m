% Angular distribution of radiation of a
% dipole interacting with a cylindrical structure

close all
clear all

load metals
lamem = 570;
nh2o = 1.33;
nal = aluminum(wavelength==lamem);
nglass = 1.52;
inner_radius = 75;
metal_thickness = 75;
dipole_pos_vec = [1e-2 5:5:70 75-1e-2];
theta = (pi/5e2:pi/2.5e2:pi)';

n0 = nh2o;
n = nh2o;
n1 = [nal nglass];
drho = 0;
nmax = 5;

for jr = 1:length(dipole_pos_vec)
    dipole_pos = dipole_pos_vec(jr);
    d0 = dipole_pos/lamem*2*pi;
    d = (inner_radius - dipole_pos)/lamem*2*pi;
    d1 = metal_thickness/lamem*2*pi;

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
    eval(['print -dpng -r300 -zbuffer adr75_' mint2str(jr) ])

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
    eval(['print -dpng -r300 -zbuffer adf75_' mint2str(jr) ])

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
    eval(['print -dpng -r300 -zbuffer adz75_' mint2str(jr) ])

    eval(['save DipoleCylADRPacBio75_' mint2str(jr) ' n0 n n1 dipole_pos_vec lamem inner_radius metal_thickness nmax theta phi adr adf adz fr ff fz'])
end

return

for j=1:36
    view([j*10 30]);
    eval(['print -dpng -r300 -zbuffer adr' mint2str(j,2)]);
end
