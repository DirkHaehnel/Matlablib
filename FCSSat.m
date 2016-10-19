if 1
    close all
    clear all
    
    sat = 10.^(-3:0.01:1);
    tfv = 3.2; % lifetime vector
    Tpulse = 0.08;
    Trep = 25;
    w = 0.3; % lateral extension of 3D-Gaussian volume
    pinhole = 50/60; % pinhole radius in object space
    lamem = 0.67;
    lamex = 0.63;
    drho = w/25; dz = drho;
    rho = 0:drho:5*w;
    z = (0.5:7*w/dz)*dz;
    [zz,rr] = meshgrid(z,rho);
    cf = D3CEF(rr/lamem*2*pi,zz/lamem*2*pi,pinhole/lamem*2*pi,1.2);
    ex = GaussianLaser(rr/lamem*2*pi,zz/lamem*2*pi,w/lamem*2*pi);
    volx = cf.*ex;    

    rhom = rho'*ones(1,length(z));
    z = ones(length(z),1)*z;
    rho = rho'*ones(1,length(rho));
    block = ones(size(rho(:,1)))*rho(:,1)';
    
    Nsub = 3;
    lam = 2^(1/Nsub);
    timewin = 1e6;
    Ncasc = ceil(log2(timewin));
    autotime = lam.^(0:Ncasc*Nsub-1)';
    dif = 5e-5; % diffusion constant in mum^2/time unit
    modres = autotime*zeros(1,length(sat));
    
    kiscv = [0 1 10 1e2]; % vector of intersystem crossing rates
    
    for jo = 1:length(tfv)
        tf = tfv(jo);
        for jk = 1:length(kiscv)
            kisc = kiscv(jk);
            for js = 1:length(sat)
                f0x = PulsedExcitation(sat(js)*volx,1/tf,Tpulse,Trep,kisc); 
                fx = f0x;
                % integrating the diffusion equation:
                for k = 1:length(autotime)
                    tmp = (erf((z-z'+dz/2)/sqrt(4*dif*autotime(k))) - erf((z-z'-dz/2)/sqrt(4*dif*autotime(k))))/2;
                    tmp = tmp + (erf((z+z'+dz/2)/sqrt(4*dif*autotime(k))) - erf((z+z'-dz/2)/sqrt(4*dif*autotime(k))))/2;        
                    fx = f0x*tmp;
                    
                    rr = rho(:,1)*rho(:,1)'/2/dif/autotime(k);
                    rpr = (rho.^2+rho'.^2)/4/dif/autotime(k);
                    rpr1 = rpr - rr;
                    
                    tmp = rpr;
                    tmp(rr<=700) = real(besselj(0, -i*rr(rr<=700)).*exp(-rpr(rr<=700)));
                    tmp(rr>700) = exp(-rpr1(rr>700))./sqrt(rr(rr>700)*2*pi);
                    tmp = block.*tmp/2/dif/autotime(k)*drho;
                    ss = sum(tmp')';
                    tmp(ss>1,:) = diag(1./ss(ss>1))*tmp(ss>1,:); 
                    fx = tmp*fx;
                    
                    modres(k,js) = real(sum(sum(rhom.*fx.*f0x)));
                    semilogx(autotime(1:k), modres(1:k,js)./modres(1,js)); axis([autotime([1 end])' 0 1]); drawnow;
                end
                modres(:,js) = modres(:,js)./modres(1,js);
            end
            
            p(1) = 1;
            for j=2:length(sat)
                p(j)=simplex('affinefit',1,0,[],[],[],autotime,modres(:,j),autotime,modres(:,1)); 
            end
            dres(:,jk,jo) = 1./p';
        end
    end
    
    save FCSSatPulsed1 dres sat tfv Tpulse Trep w drho autotime dif kiscv
    
    close all
    clear all
    
    sat = 10.^(-3:0.01:1);
    tfv = 3.2; % lifetime vector
    Tpulse = 25;
    Trep = 25;
    w = 0.3; % lateral extension of 3D-Gaussian volume
    pinhole = 50/60; % pinhole radius in onject space
    lamem = 0.67;
    lamex = 0.63;
    drho = w/25; dz = drho;
    rho = 0:drho:5*w;
    z = (0.5:7*w/dz)*dz;
    [zz,rr] = meshgrid(z,rho);
    cf = D3CEF(rr/lamem*2*pi,zz/lamem*2*pi,pinhole/lamem*2*pi,1.2);
    ex = GaussianLaser(rr/lamem*2*pi,zz/lamem*2*pi,w/lamem*2*pi);
    volx = cf.*ex;    
    
    rhom = rho'*ones(1,length(z));
    z = ones(length(z),1)*z;
    rho = rho'*ones(1,length(rho));
    block = ones(size(rho(:,1)))*rho(:,1)';

    Nsub = 3;
    lam = 2^(1/Nsub);
    timewin = 1e6;
    Ncasc = ceil(log2(timewin));
    autotime = lam.^(0:Ncasc*Nsub-1)';
    dif = 5e-5; % diffusion constant in mum^2/time unit
    modres = autotime*zeros(1,length(sat));
    
    kiscv = [0 1 10 1e2]; % vector of intersystem crossing rates
    
    for jo = 1:length(tfv)
        tf = tfv(jo);
        for jk = 1:length(kiscv)
            kisc = kiscv(jk);
            for js = 1:length(sat)
                f0x = PulsedExcitation(sat(js)*volx,1/tf,Tpulse,Trep,kisc); 
                fx = f0x;
                % integrating the diffusion equation:
                for k = 1:length(autotime)
                    tmp = (erf((z-z'+dz/2)/sqrt(4*dif*autotime(k))) - erf((z-z'-dz/2)/sqrt(4*dif*autotime(k))))/2;
                    tmp = tmp + (erf((z+z'+dz/2)/sqrt(4*dif*autotime(k))) - erf((z+z'-dz/2)/sqrt(4*dif*autotime(k))))/2;        
                    fx = f0x*tmp;
                    
                    rr = rho(:,1)*rho(:,1)'/2/dif/autotime(k);
                    rpr = (rho.^2+rho'.^2)/4/dif/autotime(k);
                    rpr1 = rpr - rr;
                    
                    tmp = rpr;
                    tmp(rr<=700) = real(besselj(0, -i*rr(rr<=700)).*exp(-rpr(rr<=700)));
                    tmp(rr>700) = exp(-rpr1(rr>700))./sqrt(rr(rr>700)*2*pi);
                    tmp = block.*tmp/2/dif/autotime(k)*drho;
                    ss = sum(tmp')';
                    tmp(ss>1,:) = diag(1./ss(ss>1))*tmp(ss>1,:); 
                    fx = tmp*fx;
                    
                    modres(k,js) = real(sum(sum(rhom.*fx.*f0x)));
                    semilogx(autotime(1:k), modres(1:k,js)./modres(1,js)); axis([autotime([1 end])' 0 1]); drawnow;
                end
                modres(:,js) = modres(:,js)./modres(1,js);
            end
            
            p(1) = 1;
            for j=2:length(sat)
                p(j)=simplex('affinefit',1,0,[],[],[],autotime,modres(:,j),autotime,modres(:,1)); 
            end
            dres(:,jk,jo) = 1./p';
        end
    end

    save FCSSatCW1 dres sat tfv Tpulse Trep w drho autotime dif kiscv

end    

% Test: 
% semilogx(autotime,modres(:,1),autotime,1./(1+4*dif/a^2*autotime)./sqrt(1+4*dif*autotime/b^2))

if 0
    load FCSSatPulsed1
    drespd=dres;
    load FCSSatCW1
    drecw=dres;
    load FCSSatData
    for j=1:size(dres,2) for k=1:size(dres,3) fpd(j,k)=simplex('affinefit',1e4,100,[],[],[],intpd1',1./pd',sat',drespd(:,j,k)); [errpd(j,k) cpd(:,j,k) zpd(:,j,k)]=affinefit(fpd(j,k),intpd1',1./pd',sat',drespd(:,j,k)); end; end
    for j=1:size(dres,2) for k=1:size(dres,3) fcw(j,k)=simplex('affinefit',1e4,100,[],[],[],intcw1',1./cw',sat',drescw(:,j,k)); [errcw(j,k) ccw(:,j,k) zcw(:,j,k)]=affinefit(fcw(j,k),intcw1',1./cw',sat',drescw(:,j,k)); end; end    
    j=2;k=2; plot(intpd2([1 2 4 5]),1./pd(1,:),'o',[0 sat*fpd(j)],cpd(1,j)+cpd(2,j)*[1; drespd(:,j)],intcw1,1./cw,'sr',[0 sat*fcw(k)],ccw(1,k)+ccw(2,k)*[1; drescw(:,k)]); axis([0 1000 0 3]);
    xlabel('excitation power [\muW]');
    ylabel('diffusion coefficient [a.u.]');
end

