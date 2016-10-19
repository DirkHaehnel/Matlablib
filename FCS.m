function [modres, autotime] = FCS(rho,z,volx1,voly1,volx2,voly2,d,velo,dif,bld)

close all
nn = size(volx1,3);
maxm = (nn-1)/2;
if nargin>4 && ~isempty(volx2)
    flag = true;
else
    flag = false;
end

if nargin<7 || isempty(d)
    d = 0;
end

if size(z,2)==1
    z = [];
else
    if size(rho,1)==1 || size(rho,2)==1
        rho = rho(:)*ones(1,numel(z));
        z = ones(numel(rho),1)*z(:)';
    end
end

if nargin<4 || isempty(voly1)
    voly1 = volx1;
end

Nsub = 3;
lam = 2^(1/Nsub);
timewin = 1e7;
Ncasc = ceil(log2(timewin));
autotime = lam.^(0:Ncasc*Nsub-1)';
if nargin<9 || isempty(dif)
    dif = 5e-5; % diffusion constant in mum^2/time unit
end

if nargin<10 || isempty(bld)
    bld = 1;
end

modres = autotime;
rhom = repmat(rho,[1 1 nn]);
if ~isempty(z)
    dz = z(1,2)-z(1,1); 
    z = ones(size(z,2),1)*z(1,:);
end
drho = rho(2,1)-rho(1,1); 
rho = rho(:,1)*ones(1,size(rho,1));
block = ones(size(rho(:,1)))*rho(:,1)';

f1 = volx1+voly1;
ff1 = f1;
if flag
    f2 = volx2+voly2;
    ff2 = f2;
end

if nargin>7 && ~isempty(velo)
    if length(velo)==1
        velo(2) = 0;
    end
    theta = -atan(velo(2)/velo(1));
    velo = sqrt(sum(velo.^2));
end

% integrating the diffusion equation:
for k = 1:length(autotime)
    if ~isempty(z)
        s4dt = sqrt(4*dif*autotime(k));
        tmp = (erf((z-z'+dz/2)/s4dt) - erf((z-z'-dz/2)/s4dt))/2;
        if abs(z(1,1))<=dz
            tmp = tmp + (erf((z+z'+dz/2)/s4dt) - erf((z+z'-dz/2)/s4dt))/2;
        end
        if d>0
            tmp = tmp + (erf((z+2*d-z'+dz/2)/s4dt) - erf((z+2*d-z'-dz/2)/s4dt))/2;
        end
        for k2 = 1:2*maxm+1
            ff1(:,:,k2) = f1(:,:,k2)*tmp;
            if flag	
                ff2(:,:,k2) = f2(:,:,k2)*tmp; 
            end
        end
    else
       ff1 = f1; 
       if flag
           ff2 = f2; 
       end
    end

    rr = rho(:,1)*rho(:,1)'/2/dif/autotime(k);
    rpr = (rho-rho').^2/4/dif/autotime(k);

    for k2 = 1:maxm+1
        tmp = besseli(k2-1, rr, 1).*exp(-rpr);
        tmp = block.*tmp/2/dif/autotime(k)*drho;
        ss = sum(tmp,2);
        tmp(ss>1,:) = diag(1./ss(ss>1))*tmp(ss>1,:);
        ff1(:,:,k2) = tmp*ff1(:,:,k2);
        if flag 
            ff2(:,:,k2) = tmp*ff2(:,:,k2); 
        end
        if k2>1
            ff1(:,:,maxm+k2) = tmp*ff1(:,:,maxm+k2);
            if flag 
                ff2(:,:,maxm+k2) = tmp*ff2(:,:,maxm+k2); 
            end
        end
    end
    
    if nargin>7 && ~isempty(velo)
        v1 = f1;
        v2 = f2;
        for jj=0:4*maxm
            psi = jj/(4*maxm+1)*2*pi;
            vr = velo*autotime(k);
            rr = sqrt(rho(:,1).^2 + vr^2 - 2*rho(:,1)*vr*cos(psi-theta));
            phi = pi + atan((vr*sin(theta)-rho(:,1)*sin(psi))./(vr*cos(theta)-rho(:,1)*cos(psi)));

            tmp = zeros(size(rho,1),size(z,2));
            for jr = 1:length(rr)
                if rr(jr)<max(rho(:,1))
                    tmp(jr,:) = interp1(rho(:,1),f1(:,:,1),rr(jr),'cubic');
                    for kk = 1:maxm 
                        tmp(jr,:) = tmp(jr,:) + interp1(rho(:,1),f1(:,:,kk+1),rr(jr),'cubic')*cos(kk*phi(jr)) + ... 
                            interp1(rho(:,1),f1(:,:,maxm+1+kk),rr(jr),'cubic')*sin(kk*phi(jr));
                    end
                end
            end
            v1(:,:,1) = v1(:,:,1) + tmp;
            for kk = 1:maxm
                v1(:,:,kk+1) = v1(:,:,kk+1) + 2*tmp*cos(kk*psi);
                v1(:,:,maxm+kk+1) = v1(:,:,maxm+kk+1) + 2*tmp*sin(kk*psi);
            end

            tmp = zeros(size(rho,1),size(z,2));
            for jr = 1:length(rr)
                if rr(jr)<max(rho(:,1))
                    tmp(jr,:) = interp1(rho(:,1),f2(:,:,1),rr(jr),'cubic');
                    for kk = 1:maxm 
                        tmp(jr,:) = tmp(jr,:) + interp1(rho(:,1),f2(:,:,kk+1),rr(jr),'cubic')*cos(kk*phi(jr)) + ... 
                            interp1(rho(:,1),f2(:,:,maxm+1+kk),rr(jr),'cubic')*sin(kk*phi(jr));
                    end
                end
            end
            v2(:,:,1) = v2(:,:,1) + tmp;
            for kk = 1:maxm
                v2(:,:,kk+1) = v2(:,:,kk+1) + 2*tmp*cos(kk*psi);
                v2(:,:,maxm+kk+1) = v2(:,:,maxm+kk+1) + 2*tmp*sin(kk*psi);
            end
            disp(jj)
        end
        save test
        
        if ~isempty(z)
            modres(k,1) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*ff1.*v1))));
            if flag
                modres(k,2) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*ff2.*v2))));
                modres(k,3) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*(ff1.*v2)))));
                modres(k,4) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*(ff2.*v1)))));
            end
        else
            modres(k,1) = real([2 ones(1,2*maxm)]*squeeze(sum(rhom.*ff1.*v1)));
            if flag
                modres(k,2) = real([2 ones(1,2*maxm)]*squeeze(sum(rhom.*ff2.*v2)));
                modres(k,3) = real([2 ones(1,2*maxm)]*squeeze(sum(rhom.*(ff1.*v2))));
                modres(k,4) = real([2 ones(1,2*maxm)]*squeeze(sum(rhom.*(ff2.*v1))));                
            end
        end
        if bld
            semilogx(autotime(1:k), modres(1:k,:)/modres(1)); axis([autotime([1 end])' 0 1]); drawnow;
        end

    else

        if ~isempty(z)
            modres(k,1) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*ff1.*f1))));
            if flag
                modres(k,2) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*ff2.*f2))));
                modres(k,3) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*(ff1.*f2+ff2.*f1)))))/2;

            end
        else
            modres(k,1) = real([2 ones(1,2*maxm)]*squeeze(sum(rhom.*ff1.*f1)));
            if flag
                modres(k,2) = real([2 ones(1,2*maxm)]*squeeze(sum(rhom.*ff2.*f2)));
                modres(k,3) = real([2 ones(1,2*maxm)]*squeeze(sum(rhom.*(ff1.*f2+ff2.*f1))))/2;
            end
        end
        if bld
            semilogx(autotime(1:k), modres(1:k,:)/modres(1)); axis([autotime([1 end])' 0 1]); drawnow;
        end

    end
end
modres = modres/modres(1,1);