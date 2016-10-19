function [modres, autotime] = FCS22(rho,z,volx1,voly1,volx2,voly2,d)

close all
nn = size(volx1,3);
maxm = (nn-1)/2;
if nargin<7 || isempty(d)
    d = 0;
end

if size(z,2)==1
    z = [];
end

Nsub = 3;
lam = 2^(1/Nsub);
timewin = 1e7;
Ncasc = ceil(log2(timewin));
autotime = lam.^(0:Ncasc*Nsub-1)';
dif = 5e-5; % diffusion constant in mum^2/time unit

modres = zeros(length(autotime),3);
rhom = repmat(rho,[1 1 nn]);
if ~isempty(z)
    dz = z(1,2)-z(1,1); 
    z = ones(size(z,2),1)*z(1,:);
end
drho = rho(2,1)-rho(1,1); 
rho = rho(:,1)*ones(1,size(rho,1));
block = ones(size(rho(:,1)))*rho(:,1)';

f1 = volx1+voly1;
f2 = volx2+voly2;

% integrating the diffusion equation:
for k = 1:length(autotime)
    ff1 = f1;
    ff2 = f2;
    if ~isempty(z)
        s4dt = sqrt(4*dif*autotime(k));
        tmpz = (erf((z-z'+dz/2)/s4dt) - erf((z-z'-dz/2)/s4dt))/2;
        if abs(z(1,1))<=dz
            tmpz = tmpz + (erf((z+z'+dz/2)/s4dt) - erf((z+z'-dz/2)/s4dt))/2;
        end
        if d>0
            tmpz = tmpz + (erf((z+2*d-z'+dz/2)/s4dt) - erf((z+2*d-z'-dz/2)/s4dt))/2;
        end
    end
    rr = rho(:,1)*rho(:,1)'/2/dif/autotime(k);
    rpr = (rho-rho').^2/4/dif/autotime(k);
    for k2 = 1:maxm+1
        tmpr(:,:,k2) = besseli(k2-1, rr, 1).*exp(-rpr);
        tmpr(:,:,k2) = block.*tmpr(:,:,k2)/2/dif/autotime(k)*drho;
        ss = sum(tmpr(:,:,k2),2);
        tmpr(ss>1,:,k2) = diag(1./ss(ss>1))*tmpr(ss>1,:,k2);
    end
    
    % 1. & 2. photon
    if ~isempty(z)
        for k2 = 1:2*maxm+1
            ff1(:,:,k2) = ff1(:,:,k2)*tmpz;
            ff2(:,:,k2) = ff2(:,:,k2)*tmpz;
        end
    end
    for k2 = 1:maxm+1
        ff1(:,:,k2) = tmpr(:,:,k2)*ff1(:,:,k2);
        ff2(:,:,k2) = tmpr(:,:,k2)*ff2(:,:,k2);
        if k2>1
            ff1(:,:,maxm+k2) = tmpr(:,:,k2)*ff1(:,:,maxm+k2);
            ff2(:,:,maxm+k2) = tmpr(:,:,k2)*ff2(:,:,maxm+k2);
        end
    end
    ff1 = FieldConv(ff1,f1);
    ff2 = FieldConv(ff2,f2);
    
    % 3. photon
    if ~isempty(z)
        for k2 = 1:2*maxm+1
            ff1(:,:,k2) = ff1(:,:,k2)*tmpz;
            ff2(:,:,k2) = ff2(:,:,k2)*tmpz;
        end
    end
    for k2 = 1:maxm+1
        ff1(:,:,k2) = tmpr(:,:,k2)*ff1(:,:,k2);
        ff2(:,:,k2) = tmpr(:,:,k2)*ff2(:,:,k2);
        if k2>1
            ff1(:,:,maxm+k2) = tmpr(:,:,k2)*ff1(:,:,maxm+k2);
            ff2(:,:,maxm+k2) = tmpr(:,:,k2)*ff2(:,:,maxm+k2);
        end
    end
    fc1 = FieldConv(ff1,f2);
    fc2 = FieldConv(ff2,f1);
    ff1 = FieldConv(ff1,f1);
    ff2 = FieldConv(ff2,f2);
    
    % 4. photon
    if ~isempty(z)
        for k2 = 1:2*maxm+1
            ff1(:,:,k2) = ff1(:,:,k2)*tmpz;
            ff2(:,:,k2) = ff2(:,:,k2)*tmpz;
            fc1(:,:,k2) = fc1(:,:,k2)*tmpz;
            fc2(:,:,k2) = fc2(:,:,k2)*tmpz;
        end
    end
    for k2 = 1:maxm+1
        ff1(:,:,k2) = tmpr(:,:,k2)*ff1(:,:,k2);
        ff2(:,:,k2) = tmpr(:,:,k2)*ff2(:,:,k2);
        fc1(:,:,k2) = tmpr(:,:,k2)*fc1(:,:,k2);
        fc2(:,:,k2) = tmpr(:,:,k2)*fc2(:,:,k2);
        if k2>1
            ff1(:,:,maxm+k2) = tmpr(:,:,k2)*ff1(:,:,maxm+k2);
            ff2(:,:,maxm+k2) = tmpr(:,:,k2)*ff2(:,:,maxm+k2);
            fc1(:,:,maxm+k2) = tmpr(:,:,k2)*fc1(:,:,maxm+k2);
            fc2(:,:,maxm+k2) = tmpr(:,:,k2)*fc2(:,:,maxm+k2);
        end
    end
    ff1 = ff1.*f1;
    ff2 = ff2.*f2;
    fc1 = fc1.*f2;
    fc2 = fc2.*f1;

    if ~isempty(z)
        modres(k,1) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*ff1))));
        modres(k,2) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*ff2))));
        modres(k,3) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*(fc1+fc2)))))/2;
    else
        modres(k,1) = real([2 ones(1,2*maxm)]*squeeze(sum(rhom.*ff1)));
        modres(k,2) = real([2 ones(1,2*maxm)]*squeeze(sum(rhom.*ff2)));
        modres(k,3) = real([2 ones(1,2*maxm)]*squeeze(sum(rhom.*(fc1+fc2))))/2;
    end
    semilogx(autotime(1:k), modres(1:k,:)/modres(1)); axis([autotime([1 end])' 0 1]); drawnow;
end
modres = modres/modres(1,1);



