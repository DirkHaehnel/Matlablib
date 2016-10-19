close all
n = 1:70;
theta = (0.5:300)'/300*pi; phi = (0:pi/79:2*pi); r = 20*pi; rr=20*pi*(0.5:100)/100;

if 0
    ephi = exp(1i*phi); ephi1 = 1./ephi;
    fr = zeros(length(theta), length(phi));
    ft = fr; ff = fr;
    for j=1:length(n)
        [ftm,ffm] = VSHM(j,1,r,theta,'j');
        [frn,ftn,ffn] = VSHN(j,1,r,theta,'j');
        tmp = 2*pi*1i.^(j+1)*sqrt((2*j+1)/4/pi/j/(j+1));
        fr = fr + tmp*frn*ephi;
        ft = ft + tmp*(ftm+ftn)*ephi;
        ff = ff + tmp*(ffm+ffn)*ephi;    
        [ftm,ffm] = VSHM(j,-1,r,theta,'j');
        [frn,ftn,ffn] = VSHN(j,-1,r,theta,'j');
        fr = fr + tmp*frn*ephi1;
        ft = ft + tmp*(-ftm+ftn)*ephi1;
        ff = ff + tmp*(-ffm+ffn)*ephi1;    
    end
    Rotor3D(real((sin(theta)*cos(phi)).*fr + (cos(theta)*cos(phi)).*ft - (ones(size(theta))*sin(phi)).*ff),[],1)
    % Rotor3D(real((sin(theta)*sin(phi)).*fr + (cos(theta)*sin(phi)).*ft + (ones(size(theta))*cos(phi)).*ff))
    % Rotor3D(real((cos(theta)*ones(size(phi))).*fr - (sin(theta)*ones(size(phi))).*ft))
end

if 0
    ephi = exp(i*phi); ephi1 = 1./ephi;
    fr = zeros(length(theta), length(phi));
    ft = fr; ff = fr;
    for j=1:length(n)
        [ftm,ffm] = VSHM(j,1,r,theta,'j');
        [frn,ftn,ffn] = VSHN(j,1,r,theta,'j');
        tmp = 2*pi*i.^j*sqrt((2*j+1)/4/pi/j/(j+1));
        fr = fr + tmp*frn*ephi;
        ft = ft + tmp*(ftm+ftn)*ephi;
        ff = ff + tmp*(ffm+ffn)*ephi;    
        [ftm,ffm] = VSHM(j,-1,r,theta,'j');
        [frn,ftn,ffn] = VSHN(j,-1,r,theta,'j');
        fr = fr - tmp*frn*ephi1;
        ft = ft + tmp*(ftm-ftn)*ephi1;
        ff = ff + tmp*(ffm-ffn)*ephi1;    
    end
    Rotor3D(real((sin(theta)*sin(phi)).*fr + (cos(theta)*sin(phi)).*ft + (ones(size(theta))*cos(phi)).*ff),[],1)
    % Rotor3D(real((sin(theta)*sin(phi)).*fr + (cos(theta)*sin(phi)).*ft + (ones(size(theta))*cos(phi)).*ff))
    % Rotor3D(real((cos(theta)*ones(size(phi))).*fr - (sin(theta)*ones(size(phi))).*ft))
end


if 0
    fr = zeros(length(theta), length(rr));
    ft = fr; ff = fr;
    fr1 = fr; ft1 = ft; ff1 = ff;
    theta0 = pi/4; x = cos(theta0); sx = sqrt(1-x^2);
    % theta0 = acos(0.2*i); x = cos(theta0);
    phi0 = 0;
    phi = 0;
    ptheta = 0;
    pphi = 1;
    [pm0,pd0] = mLegendre(max(n),x);
    [pm,pd] = mLegendre(max(n),cos(theta));
    for j = 1:length(n)
        for k = -j:j
            [ftm,ffm] = VSHM(j,k,rr,theta,'j',pm,pd);
            [frn,ftn,ffn] = VSHN(j,k,rr,theta,'j',pm,pd);

            tmp = 4*pi/1i.^(j+1)*sqrt((2*j+1)/4/pi/prod((j-abs(k)+1):(j+abs(k))))/j/(j+1)*exp(1i*k*(phi-phi0));
            %tmpm = tmp*(ptheta*sLegendre(j,k,x)-1i*pphi*dLegendre(j,k,x));
            %tmpn = tmp*(1i*pphi*sLegendre(j,k,x)-ptheta*dLegendre(j,k,x));
            tmpm = tmp*(ptheta*k*pm0(abs(k)+1,j+1)/sx+1i*pphi*pd0(abs(k)+1,j+1)*sx);
            tmpn = tmp*(1i*pphi*k*pm0(abs(k)+1,j+1)/sx+ptheta*pd0(abs(k)+1,j+1)*sx);            
            
            fr = fr + tmpn*frn;
            ft = ft + tmpm*ftm + tmpn*ftn;
            ff = ff + tmpm*ffm + tmpn*ffn;    
            fr1 = fr1 + (-1)^k*tmpn*frn;
            ft1 = ft1 + (-1)^k*(tmpm*ftm + tmpn*ftn);
            ff1 = ff1 + (-1)^k*(tmpm*ffm + tmpn*ffn);    
        end
        j
    end
    %pcolor(sin(theta)*rr,cos(theta)*rr,real((sin(theta)*ones(size(rr))).*fr + (cos(theta)*ones(size(rr))).*ft));
    fx = (sin(theta)*ones(size(rr))).*fr*cos(phi) + (cos(theta)*ones(size(rr))).*ft*cos(phi) - ff*sin(phi) ;
    fy = (sin(theta)*ones(size(rr))).*fr*sin(phi) + (cos(theta)*ones(size(rr))).*ft*sin(phi) + ff*cos(phi);
    fz = (cos(theta)*ones(size(rr))).*fr - (sin(theta)*ones(size(rr))).*ft;
    fx1 = (sin(theta)*ones(size(rr))).*fr1*cos(phi+pi) + (cos(theta)*ones(size(rr))).*ft1*cos(phi+pi) - ff1*sin(phi+pi) ;
    fy1 = (sin(theta)*ones(size(rr))).*fr1*sin(phi+pi) + (cos(theta)*ones(size(rr))).*ft1*sin(phi+pi) + ff1*cos(phi+pi);
    fz1 = (cos(theta)*ones(size(rr))).*fr1 - (sin(theta)*ones(size(rr))).*ft1;

    ftot = fx*(ptheta*cos(theta0)*cos(phi0) - pphi*sin(phi0)) + ...
        fy*(ptheta*cos(theta0)*sin(phi0) + pphi*cos(phi0)) - ...
        fz*ptheta*sin(theta0);
    ftot1 = fx1*(ptheta*cos(theta0)*cos(phi0) - pphi*sin(phi0)) + ...
        fy1*(ptheta*cos(theta0)*sin(phi0) + pphi*cos(phi0)) - ...
        fz1*ptheta*sin(theta0);
    
    pcolor([-sin(theta); sin(theta)]*rr,[cos(theta); cos(theta)]*rr,[real(ftot1); real(ftot)]);
    shading interp; axis image
    
    figure
    pcolor([-sin(theta); sin(theta)]*rr,[cos(theta); cos(theta)]*rr,...
        [cos(([sin(theta)*cos(pi+phi) sin(theta)*sin(pi+phi) cos(theta)]*[sin(theta0)*cos(phi0) sin(theta0)*sin(phi0) cos(theta0)]')*rr); ...
        cos(([sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]*[sin(theta0)*cos(phi0) sin(theta0)*sin(phi0) cos(theta0)]')*rr)]); 
    shading interp; axis image
end

if 0
    x = cos(pi-pi/3);
    fr = zeros(length(theta), length(rr));
    ft = fr; ff = fr;
    for j=1:length(n)
        [ftm,ffm] = VSHM(j,1,rr,theta);
        [frn,ftn,ffn] = VSHN(j,1,rr,theta);
        tmp = 2*pi/i.^(j+1)*sqrt((2*j+1)/4/pi/j/(j+1))/j/(j+1);
        tmp = -2*tmp*(sLegendre(j,1,x)+dLegendre(j,1,x));
        fr = fr + tmp*frn;
        ft = ft + tmp*(ftm+ftn);
        ff = ff + tmp*(ffm+ffn);    
%         [ftm,ffm] = VSHM(j,-1,rr,theta);
%         [frn,ftn,ffn] = VSHN(j,-1,rr,theta);
%         fr = fr + tmp*frn;
%         ft = ft + tmp*(-ftm+ftn);
%         ff = ff + tmp*(-ffm+ffn);    
    end
    pcolor(sin(theta)*rr,cos(theta)*rr,real((sin(theta)*ones(size(rr))).*fr + (cos(theta)*ones(size(rr))).*ft));
    shading interp; axis image
end

if 1
    NA = 1.2;
    thetamax = asin(1.2/1.33);
    x = cos((0.5:200)/200*thetamax);
    dx = thetamax/200;
    fr = zeros(length(theta), length(rr));
    ft = fr; fr1=fr; ft1=fr; 
    % ff = fr; ff1=fr;
    st = 0;
    rho0 = 10*pi; phi0 = 0; z0 = -10*pi;
    fac = exp(-0.0*(1-x.^2)-1i*z0*x).*sqrt(1-x.^2)*dx;
    tst2 = [];
    nm = n(end)+1;
    [pm,pd] = mLegendre(max(n),x);        
    [mm,zetam,zetam] = ndgrid(0:max(n),0:max(n),x);
    sL1 = mm.*pm./sqrt(1-zetam.^2);
    dL1 = sqrt(1-zetam.^2).*pd;
    zeta = cos(theta);
    [pm,pd] = mLegendre(max(n),zeta);        
%     [mm,zetam,zetam] = ndgrid(0:max(n),0:max(n),zeta);
%     sL = mm.*pm./sqrt(1-zetam.^2);
%     dL = sqrt(1-zetam.^2).*pd;
    for j=1:length(n)
        for k=-j:j
            [ftm,ffm] = VSHM(j,k,rr,theta,'j',pm,pd);
            [frn,ftn,ffn] = VSHN(j,k,rr,theta,'j',pm,pd);
            tmp = (2*pi)^2*1i.^(j-k)*sqrt((2*j+1)/4/pi/prod((j-abs(k)+1):(j+abs(k))))/j/(j+1);
            a = besselj(k-1,rho0*sqrt(1-x.^2)).*exp(-1i*(k-1)*phi0);   
            b = besselj(k+1,rho0*sqrt(1-x.^2)).*exp(-1i*(k+1)*phi0);
            tmpa = tmp*sum(fac.*a.*squeeze(sign(k)*sL1(abs(k)+1,j+1,:)-dL1(abs(k)+1,j+1,:))');
            tmpb = tmp*sum(fac.*b.*squeeze(sign(k)*sL1(abs(k)+1,j+1,:)+dL1(abs(k)+1,j+1,:))');            
            fr = fr + (tmpa+tmpb)*frn;
            ft = ft + (tmpa-tmpb)*ftm + (tmpa+tmpb)*ftn;
            fr1 = fr1 + (tmpa+tmpb)*frn*exp(1i*k*pi);
            ft1 = ft1 + ((tmpa-tmpb)*ftm + (tmpa+tmpb)*ftn)*exp(1i*k*pi);
        %   ff = ff + (tmpm-tmpn)*ffm + (tmpm+tmpn)*ffn;    
        end
        %pcolor([-sin(theta); sin(theta)]*rr,[cos(theta); cos(theta)]*rr,[-real((sin(theta)*ones(size(rr))).*fr1 + (cos(theta)*ones(size(rr))).*ft1); real((sin(theta)*ones(size(rr))).*fr + (cos(theta)*ones(size(rr))).*ft)]); 
        %colormap copper; shading flat; axis image; axis off; title(int2str(j))
        pcolor([-sin(theta); sin(theta)]*rr,[cos(theta); cos(theta)]*rr,[abs(fr1).^2+abs(ft1).^2; abs(fr).^2+abs(ft).^2]); 
        %pcolor([-sin(theta); sin(theta)]*rr,[cos(theta); cos(theta)]*rr,[real(fr1).^2+real(ft1).^2; real(fr).^2+real(ft).^2]); 
        colormap hot; shading flat; axis image; axis off; title(int2str(j))
        drawnow    
    end
    (pi*(3/2-x(end)-x(end)^2/2))^2 % max value
end
