figure(gcf);
NA = 1.2;
wd = 3e3;
tubelens = 180e3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = [];
d1 = [];
lamex = 0.635;
w0 = [500.9517 841.7505];
z0 = 1e7*[1.1879 2.0499];
over = [wd*NA 0 4.9e3];
ww = 4.9e3;
focpos = 10;
lamem = 0.67;
mag = tubelens/wd; %60;
av = 100/2; %pinhole radius!!!
zpin = 0; 
atf = [1.51 0];
kappa = 0; 
lt = [];
sat = 0; %0:0.1:3;
tsh = 1/exp(2);

if 0
    j=1;
    zeta = (j-1)/10*2*pi*ww^2/NA^2/wd^2;
    w0 = ww/sqrt(1+zeta^2);
    z0 = pi*w0^2*zeta/lamex;
    over = [wd*NA z0 -z0 w0 w0];
    
    [volx, voly, rho, z] = MDF([0 2], [8 12], NA, n0, n, n1, d0, d, d1, lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat);
    %voly = MDF([0 2], [8 12], NA, n0, n, n1, d0, d, d1, lamex, over([1 3 2 5 4]), focpos, lamem, mag, av, zpin, atf, kappa, lt, sat);
    
    rr = rho(rho(:,1)<rho(end,1)/sqrt(2),1); 
    rr = [-flipud(rr); rr];
    [xx,yy,zz] = meshgrid(rr,rr,z(1,:)');
    rr = sqrt(xx.^2+yy.^2); phi = atan(yy./xx);
    
    feld = interp2(z,rho,volx(:,:,1),zz,rr,'cubic');
    for jj=2:size(volx,3) 
        feld = feld + interp2(z,rho,volx(:,:,jj),zz,rr,'cubic').*cos(2*(jj-1)*phi); 
    end
    [tx, ty, tz] = surfcv(xx,yy,zz,feld,tsh*max(feld(:)));
    fill3(tx, ty, tz, sqrt(tx.^2+ty.^2))
    set(gcf,'renderer','zbuffer');    
    axis image
    axis([-0.65 0.65 -0.65 0.65 9 11])
    box
    shading interp;
    camlight
    lighting gouraud
    times16
    eval(['print -dpng -r100 tmp' mint2str(j,2)]);
    
    for j=2:11
        zeta = (j-1)/20*2*pi*ww^2/NA^2/wd^2;
        w0 = ww/sqrt(1+zeta^2);
        z0 = pi*w0^2*zeta/lamex;
        over = [wd*NA z0 -z0 w0 w0];
        
        volx = MDF([0 2], [8 12], NA, n0, n, n1, d0, d, d1, lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat);
        %voly = MDF([0 2], [8 12], NA, n0, n, n1, d0, d, d1, lamex, over([1 3 2 5 4]), focpos, lamem, mag, av, zpin, atf, kappa, lt, sat);
        
        feld = interp2(z,rho,volx(:,:,1),zz,rr,'cubic');
        for jj=2:size(volx,3) 
            feld = feld + interp2(z,rho,volx(:,:,jj),zz,rr,'cubic').*cos(2*(jj-1)*phi); 
        end
        [tx, ty, tz] = surfcv(xx,yy,zz,feld,tsh*max(feld(:)));
        fill3(tx, ty, tz, sqrt(tx.^2+ty.^2))
        set(gcf,'renderer','zbuffer');    
        axis image
        axis([-0.65 0.65 -0.65 0.65 9 11])
        box
        shading interp;
        camlight
        lighting gouraud
        times16
        eval(['print -dpng -r100 tmp' mint2str(j,2)]);
    end
end

if 0
    figure(gcf);
    for jsat=1:length(sat)
        
        for j=1:3 ff(:,:,:,j) = interp2(z,rho,volx(:,:,j,jsat),zz,rr,'cubic'); end
        feld = ff(:,:,:,1); for j=2:3 feld = feld + ff(:,:,:,j).*cos(2*(j-1)*phi); end
        
        [tx, ty, tz] = surfcv(xx,yy,zz,feld,max(feld(:))/exp(2));
        fill3(tx, ty, tz, sat(jsat)); caxis([0 max(sat)])
        set(gcf,'renderer','zbuffer');    
        axis image
        box
        shading flat;
        camlight
        lighting gouraud
        times16
        
        %axis([-0.5 0.5 -0.5 0.5 8.5 11.5]);
        axis([-0.7 0.7 -0.7 0.7 9.3 11.5]);
        text(0.9, -0.9, 11, ['\zeta = ' mnum2str(sat(jsat),1,1)],'fontname','times','fontsize',16);
        drawnow
        %    axis vis3d
        %     h=xlabel('\rho'); set(h,'position',sign(get(h,'position')).*[0 2.5 0]); 
        %     h=ylabel('\itz'); set(h,'position',sign(get(h,'position')).*[1.5 0 0]); 
        %     set(gca,'xtick', [-0.5 0 0.5],'xticklabel',[0.5 0 0.5],'ytick',[-1 0 1]);
        
        % figure(gcf); for j=1:360/5 
        %     camorbit(-5,0); 
        %     drawnow
        %     eval(['print -dpng -r100 tmp' mint2str(j,2)]);
        %end
        eval(['print -dpng -r100 tmp' mint2str(jsat,2)]);
    end
end

if 0
    figure(gcf);
    axis vis3d        
    for j=1:360/5 
        camorbit(-5,0); 
        drawnow
        eval(['print -dpng -r100 tmp' mint2str(j,2)]);
    end
end

if 1
    j=1;
    %n = 1.3; n1 = n;    
    atf = [1.51 (j-6)];
    [volx, voly, rho, z] = MDF([0 2], [8 12], NA, n0, n, n1, d0, d, d1, lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat);
        
    rr = rho(rho(:,1)<rho(end,1)/sqrt(2),1); 
    rr = [-flipud(rr); rr];
    [xx,yy,zz] = meshgrid(rr,rr,z(1,:)');
    rr = sqrt(xx.^2+yy.^2); phi = atan(yy./xx);
    
    feld = interp2(z,rho,volx(:,:,1),zz,rr,'cubic');
    for jj=2:size(volx,3) 
        feld = feld + interp2(z,rho,volx(:,:,jj),zz,rr,'cubic').*cos(2*(jj-1)*phi); 
    end
    [tx, ty, tz] = surfcv(xx,yy,zz,feld,tsh*max(feld(:)));
    fill3(tx, ty, tz, sqrt(tx.^2+ty.^2))
    set(gcf,'renderer','zbuffer');    
    axis image    
    axis([-0.65 0.65 -0.65 0.65 8 12])
    box
    shading interp;
    camlight
    lighting gouraud
    times16
    %text(0.75,-0.75,10,['\itn\rm = ' mnum2str(n,1,2)],'FontName','Times','FontSize',16)    
    text(0.75,-0.75,10,['\delta = ' mnum2str(atf(2),2,1) ' \mum'],'FontName','Times','FontSize',16)    
    drawnow
    eval(['print -dpng -r100 csld' mint2str(j,2)]);
    
    for j=2:11
        %n = 1.3 + (j-1)/100; n1 = n;            
        atf = [1.51 (j-6)];
        volx = MDF([0 2], [8 12], NA, n0, n, n1, d0, d, d1, lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat);
                
        feld = interp2(z,rho,volx(:,:,1),zz,rr,'cubic');
        for jj=2:size(volx,3) 
            feld = feld + interp2(z,rho,volx(:,:,jj),zz,rr,'cubic').*cos(2*(jj-1)*phi); 
        end
        [tx, ty, tz] = surfcv(xx,yy,zz,feld,tsh*max(feld(:)));
        fill3(tx, ty, tz, sqrt(tx.^2+ty.^2))
        set(gcf,'renderer','zbuffer');    
        axis image        
        axis([-0.65 0.65 -0.65 0.65 8 12])
        box
        shading interp;
        camlight
        lighting gouraud
        times16
        %text(0.75,-0.75,10,['\itn\rm = ' mnum2str(n,1,2)],'FontName','Times','FontSize',16)
        text(0.75,-0.75,10,['\delta = ' mnum2str(atf(2),2,1) ' \mum'],'FontName','Times','FontSize',16)    
        drawnow
        eval(['print -dpng -r100 csld' mint2str(j,2)]);
    end
end