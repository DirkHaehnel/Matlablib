if 0
    a=1;
    rho=1:0.1:3;
    plot(rho,(rho.^2-a^2)/2/a,-rho,(rho.^2-a^2)/2/a,'r')
    axis off
    axis image
end

if 0
    theta = -pi/2:pi/500:pi/2;
    [v,pc,ps] = DipoleL(theta,0,1.52,1,1,[],0,[]);
    plot(sin(theta').*abs(v+sign(eps+theta').*pc).^2,-cos(theta').*abs(v+sign(eps+theta').*pc).^2)
    axis image
    axis off
    hold on
    plot([0 sin(asin(1/1.52))]*10,[0 -cos(asin(1/1.52))]*10)
    plot([0 -sin(asin(1/1.52))]*10,[0 -cos(asin(1/1.52))]*10)
    hold off
end

if 0
    nm = 1;
    ng = 1.52;
    theta_c = asin(nm/ng);
    theta = (theta_c:pi/4e3:pi/2)'; % polar angle of emission; only SAF is considered
    %theta = (pi/6:pi/4e3:theta_c)'; % polar angle of emission; only subcritical emission is considered
    %theta = (pi/6:pi/1e4:pi/2)'; % polar angle of emission; both subcritical and SAF
    
    z = 0;
    [v,pc,ps] = DipoleL(theta,z,ng,nm,nm,[],z,[]);
    
    phase = unwrap(angle(ps));
    grad = gradient(phase)/mean(diff(theta));
    zp = -grad./sin(theta)/ng;
    zp(1) = zp(2);
    zp(end) = zp(end-1);
    
    phase = unwrap(angle(v));
    grad = gradient(phase)/mean(diff(theta));
    zv = -grad./sin(theta)/ng;
    zv(1) = zv(2);
    zv(end) = zv(end-1);
    
    semilogy(theta/pi*180, [zp zv], theta_c/pi*180*[1 1], [1e-1 1e2])
    legend('parallel','vertical')
    xlabel('emission angle (°)');
    ylabel('\itz\rm_0(\lambda)')
    
    
end

if 0
    clear finalx finaly finalz
    rvv = (0:5)*2*pi;
    for jrv = 1:length(rvv)
        rv = rvv(jrv);
        ParaboloidMDF;
        finalx(:,:,jrv) = fx;
        finaly(:,:,jrv) = fy;
        finalz(:,:,jrv) = fz;
    end
    finalx = cat(1,finalx(end:-1:2,:,:),finalx);
    finaly = cat(1,finaly(end:-1:2,:,:),finaly);
    finalz = cat(1,finalz(end:-1:2,:,:),finalz);
    close all
    CombineImages(cat(3,finalx,finaly,finalz),3,6);
end

if 1
    clear finalx finalz
    focposv = (0:5)*8*2*pi;
    for jrv = 1:length(focposv)
        focpos = focposv(jrv);
        ParaboloidMDF;
        finalx(:,:,jrv) = fx/max(fx(:));
        finalz(:,:,jrv) = fz/max(fz(:));
    end
    finalx = cat(1,finalx(end:-1:2,:,:),finalx);
    finalz = cat(1,finalz(end:-1:2,:,:),finalz);
    close all
    CombineImages(cat(3,finalx,finalz),2,6);
end

if 0
    clear finalx 
    focpos = 250;
    ParaboloidMDF;

    thetav = linspace(0,pi/2,6);
    for jrv=1:length(thetav)
        theta = thetav(jrv);
        finalx(:,:,jrv) = -real((sin(theta)*exx+cos(theta)*exz).*conj(sin(theta)*byx+cos(theta)*byz)-(sin(theta)*eyx+cos(theta)*eyz).*conj(sin(theta)*bxx+cos(theta)*bxz));
        finalx(:,:,jrv) = finalx(:,:,jrv)/max(max(finalx(:,:,jrv)));
    end
    finalx = cat(1,finalx(end:-1:2,:,:),finalx);
    close all
    CombineImages(finalx,1,6);
end

