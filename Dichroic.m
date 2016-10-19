figure(gcf);
clear all

if 0 % calculation of the reflectivity for varying thickness of dichroic mirror layers
    for j=0:0
        dlam = 1.37 + 0.001*j;
        
        n1 = 1.5;
        n2 = 2.3;
        
        n = [n1 n2];
        d = 2*pi*dlam./n;
        
        n = repmat(n,1,100);
        d = repmat(d,1,100);
        
        phi=0:pi/5000:pi/2; 
        [rp,rs,tp,ts] = fresnel(cos(phi),[1 n n1],d);
        plot(phi/pi*180,real(rp),'r',phi/pi*180,imag(rp),'r--',phi/pi*180,-real(rs),'b',phi/pi*180,-imag(rs),'b--',[45 45],[-1 1],':k','linewidth',1)
        %plot(phi/pi*180,abs(rp).^2,'r',phi/pi*180,abs(rs).^2,'b',[45 45],[0 1],':k','linewidth',1); axis ([0 90 0 1.002])
        times16
        xlabel('angle [°]'); ylabel('reflexion coefficients');
        title(['layer thickness = ' mnum2str(dlam,1,3)]);
        legend({'Re(\itT_p\rm)','Im(\itT_p\rm)','Re(\itT_s\rm)','Im(\itT_s\rm)'},-1)
        %legend({'|\itT_p\rm|^2','|\itT_s\rm|^2'},-1)
        drawnow;
        %eval(['print -dpng tmp' mint2str(j,3)]);
    end
end

if 1 % calculation of reflectivity at 45 degree for varying wavelength
    NLayer = 200;
    lam0 = 500;
    phase = 1.3;
    phi = pi/4;
    n1 = 1.5;
    n2 = 2.3;
    n = [n1 n2];
    n = repmat(n,1,NLayer);
    lam = 300:800;
    for j=1:length(lam)
        dlam = phase*pi*lam(j)/lam0;
        d = 2*pi*dlam./[n1 n2];
        d = repmat(d,1,NLayer);
        [rp(j),rs(j)] = fresnel(cos(phi),[1 n n1],d);
    end
    plot(lam,real(rp),'r',lam,imag(rp),'r--',lam,-real(rs),'b',lam,-imag(rs),'b--','linewidth',1)
    times16
    xlabel('wavelength [nm]'); ylabel('reflexion coefficients');
    title(['layer thickness = \lambda/' mint2str(1/phase,1)]);
    legend({'Re(\itT_p\rm)','Im(\itT_p\rm)','Re(\itT_s\rm)','Im(\itT_s\rm)'},-1)
    drawnow;
    figure
    plot(lam,abs(rp).^2,'r',lam,abs(rs).^2,'b','linewidth',1); axis ([lam(1) lam(end) 0 1.002])
    times16
    xlabel('wavelength [nm]'); ylabel('reflectivity');
    title(['layer thickness = \lambda/' mint2str(1/phase,1)]);
    legend({'|\itT_p\rm|^2','|\itT_s\rm|^2'},-1)
    drawnow;
end

