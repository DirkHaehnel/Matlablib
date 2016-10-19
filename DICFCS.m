function [modres, autotime, exc, mdf] = DICFCS(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposexc, pow, lamem, mag, av, focposdet, zpin, atf, kappa, lt, pulse, sat, triplet, resolution, ring, maxm);

close all

exc = DICExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposexc, pow, atf, resolution, maxm);

Nsub = 3;
lam = 2^(1/Nsub);
timewin = 1e6;
Ncasc = ceil(log2(timewin));
autotime = lam.^(0:Ncasc*Nsub-1)';
dif = 5e-5; % diffusion constant in mum^2/timeunit unit

modres = zeros(length(autotime),3,length(sat));

rho = exc.rho; z = exc.z;
rhom = repmat(rho,[1 1 2*maxm+1]);
dz = z(1,2)-z(1,1); drho = rho(2,1)-rho(1,1);
z = ones(size(z,2),1)*z(1,:);
rho = rho(:,1)*ones(1,size(rho,1));
block = ones(size(rho(:,1)))*rho(:,1)';

for jsat = 1:length(sat)

    mdf = DICExc2MDF(exc, NA, n0, n, n1, focposdet, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat(jsat), triplet);

    f1 = mdf.volx1+mdf.voly1;
    f2 = mdf.volx2+mdf.voly2;

    ff1 = f1;
    ff2 = f2;
    % integrating the diffusion equation:
    for k = 1:length(autotime)
        s4dt = sqrt(4*dif*autotime(k));
        tmp = (erf((z-z'+dz/2)/s4dt) - erf((z-z'-dz/2)/s4dt))/2;
        tmp = tmp + (erf((z+z'+dz/2)/s4dt) - erf((z+z'-dz/2)/s4dt))/2;
        if d>0
            tmp = tmp + (erf((z+2*d-z'+dz/2)/s4dt) - erf((z+2*d-z'-dz/2)/s4dt))/2;
        end
        for k2 = 1:2*maxm+1
            ff1(:,:,k2) = f1(:,:,k2)*tmp;
            ff2(:,:,k2) = f2(:,:,k2)*tmp;
        end

        rr = rho(:,1)*rho(:,1)'/2/dif/autotime(k);
        rpr = (rho.^2+rho'.^2)/4/dif/autotime(k);
        rpr1 = rpr - rr;

        for k2 = 1:maxm+1
            tmp = rpr;
            tmp(rr<=700) = real(i^(k2-1)*besselj(k2-1, -i*rr(rr<=700)).*exp(-rpr(rr<=700)));
            tmp(rr>700) = exp(-rpr1(rr>700))./sqrt(rr(rr>700)*2*pi);
            tmp = block.*tmp/2/dif/autotime(k)*drho;
            ss = sum(tmp')';
            tmp(ss>1,:) = diag(1./ss(ss>1))*tmp(ss>1,:);
            ff1(:,:,k2) = tmp*ff1(:,:,k2);
            ff2(:,:,k2) = tmp*ff2(:,:,k2);
            if k2>1
                ff1(:,:,maxm+k2) = tmp*ff1(:,:,maxm+k2);
                ff2(:,:,maxm+k2) = tmp*ff2(:,:,maxm+k2);
            end
        end

        modres(k,1,jsat) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*ff1.*f1))));
        modres(k,2,jsat) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*ff2.*f2))));
        modres(k,3,jsat) = real([2 ones(1,2*maxm)]*squeeze(sum(sum(rhom.*(ff1.*f2+ff2.*f1)))));
        tmp = modres(1:k,:,jsat);
        for j=1:3
            tmp(:,j) = modres(1:k,j,jsat)./max(modres(1:k,j,jsat));
        end
        semilogx(autotime(1:k), tmp); axis([autotime([1 end])' 0 1]); drawnow;
    end

end
modres = squeeze(modres);


