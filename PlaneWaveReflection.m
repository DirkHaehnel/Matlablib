close all
n = 1:20;
theta = (0:pi/99:pi)'; r = 3*pi; rr=13.2314*(0.5:50)/50;

if 0
    fr = zeros(length(theta), length(rr));
    ft = fr; ff = fr;
    fr1 = fr; ft1 = ft; ff1 = ff;
    theta0 = pi-pi/6; x = cos(theta0);
    % theta0 = acos(0.2*i); x = cos(theta0);
    phi0 = 0;
    ptheta = 1;
    pphi = 0;
    for j = 1:length(n)
        for k = -j:j
            [ftm,ffm] = VSHM(j,k,rr,theta,'j');
            [frn,ftn,ffn] = VSHN(j,k,rr,theta,'j');

            tmp = 4*pi*i.^(j+1)*sqrt((2*j+1)/4/pi/prod((j-abs(k)+1):(j+abs(k))))/j/(j+1)*exp(-i*k*phi0);
            tmpm = tmp*(ptheta*sLegendre(j,k,x)-i*pphi*dLegendre(j,k,x));
            tmpn = tmp*(i*pphi*sLegendre(j,k,x)-ptheta*dLegendre(j,k,x));            
            
            fr = fr + tmpn*frn;
            ft = ft + tmpm*ftm + tmpn*ftn;
            ff = ff + tmpm*ffm + tmpn*ffn;    
            fr1 = fr1 + (-1)^k*tmpn*frn;
            ft1 = ft1 + (-1)^k*(tmpm*ftm + tmpn*ftn);
            ff1 = ff1 + (-1)^k*(tmpm*ffm + tmpn*ffn);    
        end
    end
    %pcolor(sin(theta)*rr,cos(theta)*rr,real((sin(theta)*ones(size(rr))).*fr + (cos(theta)*ones(size(rr))).*ft));
    pcolor([-sin(theta); sin(theta)]*rr,[cos(theta); cos(theta)]*rr,[-real((sin(theta)*ones(size(rr))).*fr1 + (cos(theta)*ones(size(rr))).*ft1); real((sin(theta)*ones(size(rr))).*fr + (cos(theta)*ones(size(rr))).*ft)]); shading interp; axis image
    shading interp; axis image
end

if 1
    fr = zeros(length(theta), length(rr));
    ft = fr; ff = fr;
    fr1 = fr; ft1 = ft; ff1 = ff;
    theta0 = pi-pi/6; x = cos(theta0);
    % theta0 = acos(0.2*i); x = cos(theta0);
    phi0 = 0;
    ptheta = 1;
    pphi = 0;
    for k = 1:length(n)    
        [f, g] = VSHReflection(k,max(n),0,1.5,1,1,[],0,[]);
        tmpm = []; tmpn = [];
        for j = k:length(n)
            tmp = 4*pi*i.^(j+1)*sqrt((2*j+1)/4/pi/prod((j-abs(k)+1):(j+abs(k))))/j/(j+1)*exp(-i*k*phi0);
            tmpm(j-k+1) = tmp*(ptheta*sLegendre(j,k,x)-i*pphi*dLegendre(j,k,x));
            tmpn(j-k+1) = tmp*(i*pphi*sLegendre(j,k,x)-ptheta*dLegendre(j,k,x));            
        end
        tmpm = f*tmpm.';
        tmpn = g*tmpn.';
        for j = k:length(n)
            [ftm,ffm] = VSHM(j,k,rr,theta,'j');
            [frn,ftn,ffn] = VSHN(j,k,rr,theta,'j');

            fr = fr + tmpn(j-k+1)*frn;
            ft = ft + tmpm(j-k+1)*ftm + tmpn(j-k+1)*ftn;
            ff = ff + tmpm(j-k+1)*ffm + tmpn(j-k+1)*ffn;    
            fr1 = fr1 + (-1)^k*tmpn(j-k+1)*frn;
            ft1 = ft1 + (-1)^k*(tmpm(j-k+1)*ftm + tmpn(j-k+1)*ftn);
            ff1 = ff1 + (-1)^k*(tmpm(j-k+1)*ffm + tmpn(j-k+1)*ffn);    
        end
        pcolor([-sin(theta); sin(theta)]*rr,[cos(theta); cos(theta)]*rr,...
            log(abs([-real((sin(theta)*ones(size(rr))).*fr1 + (cos(theta)*ones(size(rr))).*ft1); ...
                real((sin(theta)*ones(size(rr))).*fr + (cos(theta)*ones(size(rr))).*ft)]))); shading interp; axis image
        drawnow
    end
    %pcolor(sin(theta)*rr,cos(theta)*rr,real((sin(theta)*ones(size(rr))).*fr + (cos(theta)*ones(size(rr))).*ft));
    pcolor([-sin(theta); sin(theta)]*rr,[cos(theta); cos(theta)]*rr,[-real((sin(theta)*ones(size(rr))).*fr1 + (cos(theta)*ones(size(rr))).*ft1); real((sin(theta)*ones(size(rr))).*fr + (cos(theta)*ones(size(rr))).*ft)]); shading interp; axis image
end


% pcolor([-sin(theta); sin(theta)]*rr,[cos(theta); cos(theta)]*rr,...
%     log(abs([-real((sin(theta)*ones(size(rr))).*fr1 + (cos(theta)*ones(size(rr))).*ft1); ...
%         real((sin(theta)*ones(size(rr))).*fr + (cos(theta)*ones(size(rr))).*ft)]))); shading interp; axis image

