% ParaboloidMDF

close all
clear all

nm = 1.0;
ng = 1.52;
zm = 0;
zinit = -15;

theta = (pi/2e4:pi/1e4:pi/2)';

if 0
    [v,pc,ps] = DipoleL(theta,zm,nm,nm,nm,[],zm,[]);
    phase = unwrap(angle(v));
    grad = gradient(phase)/mean(diff(theta));
    z0 = grad./sin(theta)/nm;
    z0(1) = z0(2);
    z0(end) = z0(end-1);
    
    zz=zinit;
    for j=1:length(theta)-1
        zz(j+1) = zz(j) + (-z0(j)-zz(j))*tan(theta(j))*mean(diff(theta)) - (z0(j+1)-z0(j))*sin(theta(j))^2;
    end
    rhorho = (-z0-zz').*tan(theta);
    
    
    [v,pc,ps] = DipoleL(theta,zm,ng,nm,nm,[],zm,[]);
    
    phase = unwrap(angle(v));
    grad = gradient(phase)/mean(diff(theta));
    z0 = grad./sin(theta)/ng;
    z0(1) = z0(2);
    z0(end) = z0(end-1);
    z0(z0>100) = min(z0);
    zv=zinit;
    for j=1:length(theta)-1
        zv(j+1) = zv(j) + (-z0(j)-zv(j))*tan(theta(j))*mean(diff(theta)) - (z0(j+1)-z0(j))*sin(theta(j))^2;
    end
    rhov = (-z0-zv').*tan(theta);
    
    phase = unwrap(angle(pc));
    grad = gradient(phase)/mean(diff(theta));
    z0 = grad./sin(theta)/ng;
    z0(1) = z0(2);
    z0(end) = z0(end-1);
    z0(z0>100) = min(z0);
    zc=zinit;
    for j=1:length(theta)-1
        zc(j+1) = zc(j) + (-z0(j)-zc(j))*tan(theta(j))*mean(diff(theta)) - (z0(j+1)-z0(j))*sin(theta(j))^2;
    end
    rhoc = (-z0-zc').*tan(theta);
    
    phase = unwrap(angle(ps));
    grad = gradient(phase)/mean(diff(theta));
    z0 = grad./sin(theta)/ng;
    z0(1) = z0(2);
    z0(end) = z0(end-1);
    z0(z0>100) = min(z0);
    zs=zinit;
    for j=1:length(theta)-1
        zs(j+1) = zs(j) + (-z0(j)-zs(j))*tan(theta(j))*mean(diff(theta)) - (z0(j+1)-z0(j))*sin(theta(j))^2;
    end
    rhos = (-z0-zs').*tan(theta);
    
    
    %plot([-rhov(end:-1:1); rhov],[zv(end:-1:1) zv],[-rhoc(end:-1:1); rhoc],[zc(end:-1:1) zc],[-rhos(end:-1:1); rhos],[zs(end:-1:1) zs],[-rhorho(end:-1:1); rhorho],[zz(end:-1:1) zz])
    plot([-rhov(end:-1:1); rhov],[zv(end:-1:1) zv],[-rhorho(end:-1:1); rhorho],[zz(end:-1:1) zz])
    axis image
    axis([-1.1*max(rhov(zv<=0)) 1.1*max(rhov(zv<=0)) min(zv) 0*max(-z0(200:200:1e3))])
    hold on
    fac = 15/cos(asin(nm/ng));
    plot([0 sin(asin(nm/ng))]*fac,[0 -cos(asin(nm/ng))]*fac,'g')
    plot([0 -sin(asin(nm/ng))]*fac,[0 -cos(asin(nm/ng))]*fac,'g')
    hold off
    
    if 0
        close all
        rho = 1;
        for j=1:10:length(theta)
            line([0 rho*sin(theta(j))],[-z0(j) -z0(j)-rho*cos(theta(j))])
        end
        line([-1.1*max(rho(z<=0)) 1.1*max(rho(z<=0))],[0 0],'color','c')
    end
    
    if 0
        plot(theta(theta>asin(nm/ng))/pi*180,-z0(theta>asin(nm/ng)))
        ax = axis;
        plot(theta(theta>asin(nm/ng))/pi*180,-z0(theta>asin(nm/ng)),asin(nm/ng)/pi*180*[1 1],[ax(3) ax(4)])
    end
end

if 0
    [v,pc,ps] = DipoleL(theta,zm,ng,nm,nm,[],zm,[]);
    [rho,z] = meshgrid(pi/60:pi/30:10*pi,0:pi/30:10*pi);
    fv = 0*rho; fp = fv;
    for j=1:length(theta)
        fv = fv + sin(theta(j))*v(j)*besselj(0,ng*rho*sin(theta(j))).*exp(1i*ng*z*cos(theta(j)));
        fp = fp + 1i*sin(theta(j))*pc(j)*besselj(1,ng*rho*sin(theta(j))).*exp(1i*ng*z*cos(theta(j)));
    end
    pcolor([-fliplr(rho) rho]/2/pi,-[z z]/2/pi,angle([fliplr(fv) fv])); 
    colorbar
    shading interp; 
    axis image
    set(gca,'yticklabel',abs(str2num(get(gca,'yticklabel'))));
    xlabel('\itx\rm (\lambda)'); ylabel('\itz\rm (\lambda)');

    pcolor([-fliplr(rho) rho]/2/pi,-[z z]/2/pi,angle([fliplr(fp) fp])); 
    colorbar
    shading interp; 
    axis image
    set(gca,'yticklabel',abs(str2num(get(gca,'yticklabel'))));
    xlabel('\itx\rm (\lambda)'); ylabel('\itz\rm (\lambda)');

    pcolor([-fliplr(rho) rho]/2/pi,-[z z]/2/pi,real([fliplr(fv) fv])/max(max(abs(real(fv))))); 
    colorbar
    shading interp; 
    axis image
    set(gca,'yticklabel',abs(str2num(get(gca,'yticklabel'))));
    xlabel('\itx\rm (\lambda)'); ylabel('\itz\rm (\lambda)');

    pcolor([-fliplr(rho) rho]/2/pi,-[z z]/2/pi,real([fliplr(fp) fp])/max(max(abs(real(fp))))); 
    colorbar
    shading interp; 
    axis image
    set(gca,'yticklabel',abs(str2num(get(gca,'yticklabel'))));
    xlabel('\itx\rm (\lambda)'); ylabel('\itz\rm (\lambda)');

    if 0
        psi = pi/40:pi/20:2*pi;
        for j=1:length(psi)
            pcolor([-fliplr(rho) rho]/2/pi,-[z z]/2/pi,angle(exp(-1i*psi(j))*[fliplr(fv) fv])); 
            shading interp;
            axis image
            set(gca,'yticklabel',abs(str2num(get(gca,'yticklabel'))));
            xlabel('\itx\rm (\lambda)'); ylabel('\itz\rm (\lambda)');
            mov(j) = getframe;
        end
    end

    if 0
        psi = pi/40:pi/20:2*pi;
        for j=1:length(psi)
            pcolor([-fliplr(rho) rho]/2/pi,-[z z]/2/pi,angle(exp(-1i*psi(j))*[fliplr(fp) fp]));
            shading interp;
            axis image
            set(gca,'yticklabel',abs(str2num(get(gca,'yticklabel'))));
            xlabel('\itx\rm (\lambda)'); ylabel('\itz\rm (\lambda)');
            mov(j) = getframe;
        end
    end

end

