% Quantum interference of single-molecule fluorescence trhough near-field stimulated emission 

close all
focus = 0.1:0.1:4;

mag = 100;
NA = 1.25;
lambda = 0.57;
focus = 1;
maxr = 2;
dist = 0;
om1 = 0*pi/2; 
om2 = 0*pi/2;


if 1 % simple dipole
    [intx inty int, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
        SEPDipole([0 sqrt(maxr^2+(maxr+dist)^2)], 0, NA, 1.51, 1., 1., [], 0, [], lambda, mag, focus);

    drho = diff(rho(1:2));
    psi1 = pi/4;
    psi2 = pi/4;
    om1 = 0;
    omv = 0; %0:pi/30:2*pi-pi/30;
    [x,y] = meshgrid(-mag*maxr:drho:mag*maxr,-mag*maxr:drho:mag*maxr);
    p1 = angle(x+1i*y)-psi1;
    p2 = angle(x+1i*y)-psi2;
    r1 = sqrt(x.^2+y.^2); r2 = r1;
    for j=1:length(omv)
        om2 = omv(j);
        
        ex1 = cos(om1)*(interp1(rho,fxx0,r1,'cubic')+cos(2*p1).*interp1(rho,fxx2,r1,'cubic'))+sin(om1)*cos(p1).*interp1(rho,fxz,r1,'cubic');
        ex2 = cos(om2)*(interp1(rho,fxx0,r2,'cubic')+cos(2*p2).*interp1(rho,fxx2,r2,'cubic'))+sin(om2)*cos(p2).*interp1(rho,fxz,r2,'cubic');
        bx1 = cos(om1)*(interp1(rho,byx0,r1,'cubic')+cos(2*p1).*interp1(rho,byx2,r1,'cubic'))+sin(om1)*cos(p1).*interp1(rho,byz,r1,'cubic');
        bx2 = cos(om2)*(interp1(rho,byx0,r2,'cubic')+cos(2*p2).*interp1(rho,byx2,r2,'cubic'))+sin(om2)*cos(p2).*interp1(rho,byz,r2,'cubic');
        ey1 = cos(om1)*sin(2*p1).*interp1(rho,fxx2,r1,'cubic')+sin(om1)*sin(p1).*interp1(rho,fxz,r1,'cubic');
        ey2 = cos(om2)*sin(2*p2).*interp1(rho,fxx2,r2,'cubic')+sin(om2)*sin(p2).*interp1(rho,fxz,r2,'cubic');    
        by1 = cos(om1)*sin(2*p1).*interp1(rho,byx2,r1,'cubic')+sin(om1)*sin(p1).*interp1(rho,byz,r1,'cubic');
        by2 = cos(om2)*sin(2*p2).*interp1(rho,byx2,r2,'cubic')+sin(om2)*sin(p2).*interp1(rho,byz,r2,'cubic');
        
        interf = real((ex1*cos(psi1)-ey1*sin(psi1)+ex2*cos(psi2)-ey2*sin(psi2)).*conj(bx1*cos(psi1)-by1*sin(psi1)+bx2*cos(psi2)-by2*sin(psi2)) + ...
            (ey1*cos(psi1)+ex1*sin(psi1)+ey2*cos(psi2)+ex2*sin(psi2)).*conj(by1*cos(psi1)+bx1*sin(psi1)+by2*cos(psi2)+bx2*sin(psi2)));
        int1 = real(ex1.*conj(bx1) + ey1.*conj(by1));
        int2 = real(ex2.*conj(bx2) + ey2.*conj(by2));    

        subplot(1,3,1);
        rad = 2;
        pfeil([1 0 0] - [0 1 0]*cos(om2) + [0 0 1]*sin(om2), [1 0 0] - [0 1.5 0]*cos(om2) + [0 0 1.5]*sin(om2),[],[],[],[],rad); 
        hold on; 
        pfeil([1 0 0] - [0 1 0]*cos(om2) + [0 0 1]*sin(om2), [1 0 0] - [0 0.5 0]*cos(om2) + [0 0 0.5]*sin(om2),[],[],[],[],rad); 
        pfeil([0 0 0], [0.5 0 0],[],[],[],[],rad);
        pfeil([0 0 0], [-0.5 0 0],[],[],[],[],rad);
        hold off; 
        caxis([-1 1.5]); colormap jet
        view([20 10])
        patch([-0.5 -0.5 1.1 1.1],[-1.5 1.5 1.5 -1.5],-2.*[1 1 1 1],0.9*[1 1 1 1],'edgecolor','none','facealpha',0.3,'edgealpha',0.3)
        axis([-0.5 1.1 -1.5 1.5 -5 4])
        axis off
        axis image
        camlight
        
        subplot(1,3,2); pcolor(x,y,int1+int2); 
        axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); axis image; shading interp; axis off; 
        title('no interference')
        subplot(1,3,3); pcolor(x,y,interf); 
        axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); axis image; shading interp; axis off
        title('with interference')    
        drawnow
        colormap hot
        %eval(['print -dpng tmpd' mint2str(j,2)]);
    end
end

if 0 % simple dipole
    [intx inty int, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
        SEPDipole([0 sqrt(maxr^2+(maxr+dist)^2)], 0, NA, 1.51, 1., 1., [], 0, [], lambda, mag, focus);

    drho = diff(rho(1:2));
    psi1 = 0;
    psi2 = pi/2;
    om1 = 0;
    omv = 0:pi/30:pi-pi/30;
    [x,y] = meshgrid(-mag*maxr:drho:mag*maxr,-mag*maxr:drho:mag*maxr);
    p1 = angle(x+i*y)-psi1;
    p2 = angle(x+i*y)-psi2;
    r1 = sqrt(x.^2+y.^2); r2 = r1;
    for j=1:length(omv)
        om2 = omv(j);
        
        ex1 = cos(om1)*(interp1(rho,fxx0,r1,'cubic')+cos(2*p1).*interp1(rho,fxx2,r1,'cubic'))+sin(om1)*cos(p1).*interp1(rho,fxz,r1,'cubic');
        ex2 = cos(om2)*(interp1(rho,fxx0,r2,'cubic')+cos(2*p2).*interp1(rho,fxx2,r2,'cubic'))+sin(om2)*cos(p2).*interp1(rho,fxz,r2,'cubic');
        bx1 = cos(om1)*(interp1(rho,byx0,r1,'cubic')+cos(2*p1).*interp1(rho,byx2,r1,'cubic'))+sin(om1)*cos(p1).*interp1(rho,byz,r1,'cubic');
        bx2 = cos(om2)*(interp1(rho,byx0,r2,'cubic')+cos(2*p2).*interp1(rho,byx2,r2,'cubic'))+sin(om2)*cos(p2).*interp1(rho,byz,r2,'cubic');
        ey1 = cos(om1)*sin(2*p1).*interp1(rho,fxx2,r1,'cubic')+sin(om1)*sin(p1).*interp1(rho,fxz,r1,'cubic');
        ey2 = cos(om2)*sin(2*p2).*interp1(rho,fxx2,r2,'cubic')+sin(om2)*sin(p2).*interp1(rho,fxz,r2,'cubic');    
        by1 = cos(om1)*sin(2*p1).*interp1(rho,byx2,r1,'cubic')+sin(om1)*sin(p1).*interp1(rho,byz,r1,'cubic');
        by2 = cos(om2)*sin(2*p2).*interp1(rho,byx2,r2,'cubic')+sin(om2)*sin(p2).*interp1(rho,byz,r2,'cubic');
        
        interf = real((ex1*cos(psi1)-ey1*sin(psi1)+ex2*cos(psi2)-ey2*sin(psi2)).*conj(bx1*cos(psi1)-by1*sin(psi1)+bx2*cos(psi2)-by2*sin(psi2)) + ...
            (ey1*cos(psi1)+ex1*sin(psi1)+ey2*cos(psi2)+ex2*sin(psi2)).*conj(by1*cos(psi1)+bx1*sin(psi1)+by2*cos(psi2)+bx2*sin(psi2)));
        int1 = real(ex1.*conj(bx1) + ey1.*conj(by1));
        int2 = real(ex2.*conj(bx2) + ey2.*conj(by2));    

        subplot(2,2,1); pcolor(x,y,int1); 
        axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); axis image; shading interp; axis off; 
        title('first molecule')
        subplot(2,2,2); pcolor(x,y,int2); 
        axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); axis image; shading interp; axis off
        title('second molecule')    
        subplot(2,2,3); pcolor(x,y,int1+int2); 
        axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); axis image; shading interp; axis off; 
        title('no interference')
        subplot(2,2,4); pcolor(x,y,interf); 
        axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); axis image; shading interp; axis off
        title('with interference')    
        drawnow
        colormap hot
        eval(['print -dpng tmpd' mint2str(j,2)]);
    end
end


if 0 % simple dipole
    [intx inty int, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
        SEPDipole([0 sqrt(maxr^2+(maxr+dist)^2)], 0, NA, 1.51, 1., 1., [], 0, [], lambda, mag, focus);

    drho = diff(rho(1:2));
    psi1 = 0;
    psiv = 0:pi/30:pi*2;
    [x,y] = meshgrid(-mag*maxr:drho:mag*maxr,-mag*maxr:drho:mag*maxr);
    p1 = angle(x+i*y)-psi1;
    r1 = sqrt(x.^2+y.^2); r2 = r1;

    j=1;
    psi2 = psiv(j);
    
    p2 = angle(x+i*y)-psi2;
    
    ex1 = cos(om1)*(interp1(rho,fxx0,r1,'cubic')+cos(2*p1).*interp1(rho,fxx2,r1,'cubic'))+sin(om1)*cos(p1).*interp1(rho,fxz,r1,'cubic');
    ex2 = cos(om2)*(interp1(rho,fxx0,r2,'cubic')+cos(2*p2).*interp1(rho,fxx2,r2,'cubic'))+sin(om2)*cos(p2).*interp1(rho,fxz,r2,'cubic');
    bx1 = cos(om1)*(interp1(rho,byx0,r1,'cubic')+cos(2*p1).*interp1(rho,byx2,r1,'cubic'))+sin(om1)*cos(p1).*interp1(rho,byz,r1,'cubic');
    bx2 = cos(om2)*(interp1(rho,byx0,r2,'cubic')+cos(2*p2).*interp1(rho,byx2,r2,'cubic'))+sin(om2)*cos(p2).*interp1(rho,byz,r2,'cubic');
    ey1 = cos(om1)*sin(2*p1).*interp1(rho,fxx2,r1,'cubic')+sin(om1)*sin(p1).*interp1(rho,fxz,r1,'cubic');
    ey2 = cos(om2)*sin(2*p2).*interp1(rho,fxx2,r2,'cubic')+sin(om2)*sin(p2).*interp1(rho,fxz,r2,'cubic');    
    by1 = cos(om1)*sin(2*p1).*interp1(rho,byx2,r1,'cubic')+sin(om1)*sin(p1).*interp1(rho,byz,r1,'cubic');
    by2 = cos(om2)*sin(2*p2).*interp1(rho,byx2,r2,'cubic')+sin(om2)*sin(p2).*interp1(rho,byz,r2,'cubic');
    
    interf = real((ex1*cos(psi1)-ey1*sin(psi1)+ex2*cos(psi2)-ey2*sin(psi2)).*conj(bx1*cos(psi1)-by1*sin(psi1)+bx2*cos(psi2)-by2*sin(psi2)) + ...
        (ey1*cos(psi1)+ex1*sin(psi1)+ey2*cos(psi2)+ex2*sin(psi2)).*conj(by1*cos(psi1)+bx1*sin(psi1)+by2*cos(psi2)+bx2*sin(psi2)));
    int1 = real(ex1.*conj(bx1) + ey1.*conj(by1));
    int2 = real(ex2.*conj(bx2) + ey2.*conj(by2));    
    
    mm = max([max(int1(:))*2 max(interf(:))]);
    
    subplot(2,2,1); pcolor(x,y,int1); 
    axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); caxis([0 mm]); axis image; shading interp; axis off; 
    title('first molecule')
    subplot(2,2,2); pcolor(x,y,int2); 
    axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); caxis([0 mm]); axis image; shading interp; axis off
    title('second molecule')    
    subplot(2,2,3); pcolor(x,y,int1+int2); 
    axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); caxis([0 mm]); axis image; shading interp; axis off; 
    title('no interference')
    subplot(2,2,4); pcolor(x,y,interf); 
    axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); caxis([0 mm]); axis image; shading interp; axis off
    title('with interference')    
    drawnow
    colormap hot
    eval(['print -dpng tmpd' mint2str(j,2)]);
    for j=2:length(psiv)
        psi2 = psiv(j);
        
        p2 = angle(x+i*y)-psi2;
        
        ex1 = cos(om1)*(interp1(rho,fxx0,r1,'cubic')+cos(2*p1).*interp1(rho,fxx2,r1,'cubic'))+sin(om1)*cos(p1).*interp1(rho,fxz,r1,'cubic');
        ex2 = cos(om2)*(interp1(rho,fxx0,r2,'cubic')+cos(2*p2).*interp1(rho,fxx2,r2,'cubic'))+sin(om2)*cos(p2).*interp1(rho,fxz,r2,'cubic');
        bx1 = cos(om1)*(interp1(rho,byx0,r1,'cubic')+cos(2*p1).*interp1(rho,byx2,r1,'cubic'))+sin(om1)*cos(p1).*interp1(rho,byz,r1,'cubic');
        bx2 = cos(om2)*(interp1(rho,byx0,r2,'cubic')+cos(2*p2).*interp1(rho,byx2,r2,'cubic'))+sin(om2)*cos(p2).*interp1(rho,byz,r2,'cubic');
        ey1 = cos(om1)*sin(2*p1).*interp1(rho,fxx2,r1,'cubic')+sin(om1)*sin(p1).*interp1(rho,fxz,r1,'cubic');
        ey2 = cos(om2)*sin(2*p2).*interp1(rho,fxx2,r2,'cubic')+sin(om2)*sin(p2).*interp1(rho,fxz,r2,'cubic');    
        by1 = cos(om1)*sin(2*p1).*interp1(rho,byx2,r1,'cubic')+sin(om1)*sin(p1).*interp1(rho,byz,r1,'cubic');
        by2 = cos(om2)*sin(2*p2).*interp1(rho,byx2,r2,'cubic')+sin(om2)*sin(p2).*interp1(rho,byz,r2,'cubic');
        
        interf = real((ex1*cos(psi1)-ey1*sin(psi1)+ex2*cos(psi2)-ey2*sin(psi2)).*conj(bx1*cos(psi1)-by1*sin(psi1)+bx2*cos(psi2)-by2*sin(psi2)) + ...
            (ey1*cos(psi1)+ex1*sin(psi1)+ey2*cos(psi2)+ex2*sin(psi2)).*conj(by1*cos(psi1)+bx1*sin(psi1)+by2*cos(psi2)+bx2*sin(psi2)));
        int1 = real(ex1.*conj(bx1) + ey1.*conj(by1));
        int2 = real(ex2.*conj(bx2) + ey2.*conj(by2));    
        
        subplot(2,2,1); pcolor(x,y,int1); 
        axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); caxis([0 mm]); axis image; shading interp; axis off; 
        title('first molecule')
        subplot(2,2,2); pcolor(x,y,int2); 
        axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); caxis([0 mm]); axis image; shading interp; axis off
        title('second molecule')    
        subplot(2,2,3); pcolor(x,y,int1+int2); 
        axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); caxis([0 mm]); axis image; shading interp; axis off; 
        title('no interference')
        subplot(2,2,4); pcolor(x,y,interf); 
        axis([-mag*maxr mag*maxr -mag*maxr mag*maxr]); caxis([0 mm]); axis image; shading interp; axis off
        title('with interference')    
        drawnow
        colormap hot
        eval(['print -dpng tmpd' mint2str(j,2)]);
    end
end


dist = 3;
om1 = 0*pi/2; 
om2 = 0*pi/2;
psi1 = pi/4;
psi2 = -pi/4;
if 0 % simple dipole
    [intx inty int, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
        SEPDipole([0 sqrt(maxr^2+(maxr+dist)^2)], 0, NA, 1.51, 1., 1., [], 0, [], lambda, mag, focus);

    drho = diff(rho(1:2));
    distv = -dist:0.1:dist;
    [x,y] = meshgrid(-mag*(maxr+max(distv)/2):drho:mag*(maxr+max(distv)/2),-mag*maxr:drho:mag*maxr);
    for j=1:length(distv)
        dist = distv(j);
        
        p1 = angle((x-mag*dist/2)+i*y)-psi1;
        r1 = sqrt((x-mag*dist/2).^2+y.^2);
        p2 = angle(x+mag*dist/2+i*y)-psi2;
        r2 = sqrt((x+mag*dist/2).^2+y.^2);
        
        ex1 = cos(om1)*(interp1(rho,fxx0,r1,'cubic')+cos(2*p1).*interp1(rho,fxx2,r1,'cubic'))+sin(om1)*cos(p1).*interp1(rho,fxz,r1,'cubic');
        ex2 = cos(om2)*(interp1(rho,fxx0,r2,'cubic')+cos(2*p2).*interp1(rho,fxx2,r2,'cubic'))+sin(om2)*cos(p2).*interp1(rho,fxz,r2,'cubic');
        bx1 = cos(om1)*(interp1(rho,byx0,r1,'cubic')+cos(2*p1).*interp1(rho,byx2,r1,'cubic'))+sin(om1)*cos(p1).*interp1(rho,byz,r1,'cubic');
        bx2 = cos(om2)*(interp1(rho,byx0,r2,'cubic')+cos(2*p2).*interp1(rho,byx2,r2,'cubic'))+sin(om2)*cos(p2).*interp1(rho,byz,r2,'cubic');
        ey1 = cos(om1)*sin(2*p1).*interp1(rho,fxx2,r1,'cubic')+sin(om1)*sin(p1).*interp1(rho,fxz,r1,'cubic');
        ey2 = cos(om2)*sin(2*p2).*interp1(rho,fxx2,r2,'cubic')+sin(om2)*sin(p2).*interp1(rho,fxz,r2,'cubic');    
        by1 = cos(om1)*sin(2*p1).*interp1(rho,byx2,r1,'cubic')+sin(om1)*sin(p1).*interp1(rho,byz,r1,'cubic');
        by2 = cos(om2)*sin(2*p2).*interp1(rho,byx2,r2,'cubic')+sin(om2)*sin(p2).*interp1(rho,byz,r2,'cubic');
        
        interf = real((ex1*cos(psi1)-ey1*sin(psi1)+ex2*cos(psi2)-ey2*sin(psi2)).*conj(bx1*cos(psi1)-by1*sin(psi1)+bx2*cos(psi2)-by2*sin(psi2)) + ...
            (ey1*cos(psi1)+ex1*sin(psi1)+ey2*cos(psi2)+ex2*sin(psi2)).*conj(by1*cos(psi1)+bx1*sin(psi1)+by2*cos(psi2)+bx2*sin(psi2)));
        int1 = real(ex1.*conj(bx1) + ey1.*conj(by1));
        int2 = real(ex2.*conj(bx2) + ey2.*conj(by2));    
        
        subplot(2,2,1); pcolor(x,y,int1); 
        axis([-mag*(maxr+max(distv)/2) mag*(maxr+max(distv)/2) -mag*maxr mag*maxr]); axis image; shading interp; axis off; 
        title('first molecule')
        subplot(2,2,2); pcolor(x,y,int2); 
        axis([-mag*(maxr+max(distv)/2) mag*(maxr+max(distv)/2) -mag*maxr mag*maxr]); axis image; shading interp; axis off
        title('second molecule')    
        subplot(2,2,3); pcolor(x,y,int1+int2); 
        axis([-mag*(maxr+max(distv)/2) mag*(maxr+max(distv)/2) -mag*maxr mag*maxr]); axis image; shading interp; axis off; 
        title('no interference')
        subplot(2,2,4); pcolor(x,y,interf); 
        axis([-mag*(maxr+max(distv)/2) mag*(maxr+max(distv)/2) -mag*maxr mag*maxr]); axis image; shading interp; axis off
        title('with interference')    
        drawnow
        colormap hot
        eval(['print -djpeg tmpd' mint2str(j,2)]);
    end
end

