lambda = 500; 
n = 1.33; 
na = 1.2; 
fd = 3e3;
Theta = asin(na/n);
a = 2*pi/2/n/sin(Theta); 
b = 2*pi/n/(1 - cos(Theta));

if 0
    [xx,zz] = meshgrid(linspace(-a,a,100), linspace(-b,b,100));
    
    muv = linspace(0,pi,200);
    
    fieldx = zeros(size(xx));
    fieldz = fieldx;
    yy = 0*xx;
    for j=1:length(muv)
        mu = muv(j);
        psiv = linspace(0,2*pi,1 + round(400*sin(mu)));
        for k=1:length(psiv)
            psi = psiv(k);
            [ex, ez] = DipoleFree(xx,yy,zz,a*sin(mu)*cos(psi),a*sin(mu)*sin(psi),b*cos(mu),n);
            fieldx = fieldx + squeeze(b*sin(mu)*(cos(psi)*ex(1,:,:) + cos(3*psi)*ex(2,:,:)) + a*cos(mu)*cos(psi)*ez(2,:,:))/sqrt(b^2*sin(mu)^2+a^2*cos(mu)^2);
            fieldz = fieldz + squeeze(-b*sin(mu)*ex(3,:,:) + a*cos(mu)*ez(1,:,:))/sqrt(b^2*sin(mu)^2+a^2*cos(mu)^2);
        end
        disp(j)
    end
    
    close all
    mpcolor(xx,zz,log(abs(fieldx).^2)); shading interp
    figure
    mpcolor(xx,zz,log(abs(fieldz).^2)); shading interp
    
    figure
    tmp=log10(abs(fieldx).^2+abs(fieldz).^2);
    contourf(xx(abs(zz(:,1))<4,abs(xx(1,:))<1), zz(abs(zz(:,1))<4,abs(xx(1,:))<1), tmp(abs(zz(:,1))<4,abs(xx(1,:))<1)-min(tmp(:)),10);
    colormap jet; colorbar; axis image
end

if 1
    rho_max = fd*na;
    [rho,phi] = meshgrid(linspace(0,rho_max,100), linspace(0,2*pi,360));
    theta = asin(rho/fd/n);
    
    muv = linspace(0,pi,200);
    
    fieldx = zeros(size(xx));
    fieldz = fieldx;
    yy = 0*xx;
    for j=1:length(muv)
        mu = muv(j);
        psiv = linspace(0,2*pi,1 + round(400*sin(mu)));
        for k=1:length(psiv)
            psi = psiv(k);
            [v, pc, ps] = DipoleL(theta,z,n0,n,n1,d0,d,d1);
            [ex, ez] = DipoleFree(xx,yy,zz,a*sin(mu)*cos(psi),a*sin(mu)*sin(psi),b*cos(mu),n);
            fieldx = fieldx + squeeze(b*sin(mu)*(cos(psi)*ex(1,:,:) + cos(3*psi)*ex(2,:,:)) + a*cos(mu)*cos(psi)*ez(2,:,:))/sqrt(b^2*sin(mu)^2+a^2*cos(mu)^2);
            fieldz = fieldz + squeeze(-b*sin(mu)*ex(3,:,:) + a*cos(mu)*ez(1,:,:))/sqrt(b^2*sin(mu)^2+a^2*cos(mu)^2);
        end
        disp(j)
    end
    
end