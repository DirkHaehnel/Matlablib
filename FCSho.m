function [modres, autotime] = FCSho(rho,z,vol,d,ord)

% program for calculating higher order ACF

close all
nn = size(vol,3);
maxm = (nn-1)/2;
if nargin<4 || isempty(d)
    d = 0;
end
if nargin<5 || isempty(ord)
    ord = 2;
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

modres = zeros(length(autotime),ord-1);
rhom = repmat(rho,[1 1 nn]);
if ~isempty(z)
    dz = z(1,2)-z(1,1); 
    z = ones(size(z,2),1)*z(1,:);
end
drho = rho(2,1)-rho(1,1); 
rho = rho(:,1)*ones(1,size(rho,1));
block = ones(size(rho(:,1)))*rho(:,1)';

f = vol;

% integrating the diffusion equation:
for k = 1:length(autotime)
    ff = f;
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

    for j=1:ord-1
        if ~isempty(z)
            for k2 = 1:2*maxm+1
                ff(:,:,k2) = ff(:,:,k2)*tmpz;
            end
        end
        for k2 = 1:maxm+1
            ff(:,:,k2) = tmpr(:,:,k2)*ff(:,:,k2);
            if k2>1
                ff(:,:,maxm+k2) = tmpr(:,:,k2)*ff(:,:,maxm+k2);
            end
        end
        if ~isempty(z)
            modres(k,j) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*ff.*f))));
        else
            modres(k,j) = real([2 ones(1,2*maxm)]*squeeze(sum(rhom.*ff.*f)));
        end
        if j<ord-1
            ff = FieldConv(ff,f);
        end
    end    
    
    semilogx(autotime(1:k), modres(1:k,:)./(ones(k,1)*modres(1,:))); axis([autotime([1 end])' 0 1]); drawnow;
end
modres = modres./modres(1,1);

