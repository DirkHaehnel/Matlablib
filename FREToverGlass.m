if 1 % FRET over glass surface in water (neglecting donor lifetime)
    close all
    lambda = 550;
    rad = 2*pi/lambda*0.1*exp(log(10/0.1)*(0:0.01:1));
    phi = (0:90)/180*pi;
    zd = 0.5*2*pi/lambda;
    n0 = 1.52;
    nd = 1.33;
    ni = 1.33;
    na = 1.33;
    n1 = 1.33;
    d0 = [];
    dd = zd;
    di = 0*pi/lambda;
    da = 0;
    d1 = [];
    
    ex = zeros(3,length(rad),length(phi));
    ex0 = ex;
    ez = zeros(2,length(rad),length(phi));
    ez0 = ez;
    
    for jr = 1:length(rad)
        for jp = 1:length(phi)
            [ex(:,jr,jp),ez(:,jr,jp),ex0(:,jr,jp),ez0(:,jr,jp)] = DipoleField(zd,rad(jr)*sin(phi(jp)),rad(jr)*cos(phi(jp)),n0,nd,ni,na,n1,d0,dd,di,da,d1);
        end
        jr
    end
    
    extotal = squeeze(abs(ex0(1,:,:) + ex(1,:,:)).^2+0.5*abs(ex0(2,:,:) + ex(2,:,:)).^2+0.5*abs(ex0(3,:,:) + ex(3,:,:)).^2);
    ex0total = squeeze(abs(ex0(1,:,:)).^2+0.5*abs(ex0(2,:,:)).^2+0.5*abs(ex0(3,:,:)).^2);
    eztotal = squeeze(abs(ez0(1,:,:) + ez(1,:,:)).^2+abs(ez0(2,:,:) + ez(2,:,:)).^2);
    ez0total = squeeze(abs(ez0(1,:,:)).^2+abs(ez0(2,:,:)).^2);
    
    cameratoolbar; surf(rad'*sin(phi)/2/pi*lambda,rad'*cos(phi)/2/pi*lambda,(2*extotal+eztotal)./(2*ex0total+ez0total)); set(gca,'linewidth',1);
    xlabel('lateral distance (nm)'); ylabel('distance from surface (nm)'); zlabel('enhancement factor  \it{Q}');

    tmp = contour(rad'*sin(phi)/2/pi*lambda,rad'*cos(phi)/2/pi*lambda,(2*extotal+eztotal)./(2*ex0total+ez0total),[1 1]);
    contourf(rad'*sin(phi)/2/pi*lambda,rad'*cos(phi)/2/pi*lambda,(2*extotal+eztotal)./(2*ex0total+ez0total)); set(gca,'linewidth',1);
    axis image
    caxis([0.85 1.1])
    %shading interp
    line(tmp(1,2:end),tmp(2,2:end),'color','w','linewidth',2)
    xlabel('lateral distance (nm)'); ylabel('distance from surface (nm)'); 
    h = colorbar
    set(h,'linewidth',1)

    
    figure
    cameratoolbar; surf(rad'*sin(phi)/2/pi*lambda,rad'*cos(phi)/2/pi*lambda,extotal./ex0total); camlight; shading interp; view(115,16); alpha(0.8)
    figure
    cameratoolbar; surf(rad'*sin(phi)/2/pi*lambda,rad'*cos(phi)/2/pi*lambda,eztotal./ez0total); camlight; shading interp; view(115,16); alpha(0.8)
    %plot(rhoa/2/pi*lambda,ez0total,rhoa/2/pi*lambda,eztotal,'o-')
    %plot(rhoa/2/pi*lambda,ex0total,rhoa/2/pi*lambda,extotal,'o-')
end

if 0
    z = (0:0.1:100)/550*2*pi; 
    [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu] = LifetimeL(z,1.52,1.33,1.33,[],0,[]);   
    z = z/2/pi*550;
    qv = 1./(4*1.33/3./(qvd+qvu));
    qp = 1./(4*1.33/3./(qpd+qpu));
    qt = 1./(4*1.33./(2*qpd+2*qpu+qvu+qvd));
    
    close
    plot(z,qv,z,qp,z,qt);
    patch([z(1:end-1);z(2:end);z(2:end);z(1:end-1)],[qp(1:end-1); qp(1:end-1); qv(1:end-1); qv(1:end-1)],'g','EdgeColor','none','FaceAlpha',0.3)
    xlabel('distance (nm)')
    ylabel('radiation enhancement    \gamma')
    axis([-0.001 100.0001 1 1.4])
    set(gca,'LineWidth',1)
    legend({'  \perp', '  ||', '\langle \cdot \rangle_\Omega'})
    print -djpeg -r600 EmissionEnhancementOverSurface
    
end