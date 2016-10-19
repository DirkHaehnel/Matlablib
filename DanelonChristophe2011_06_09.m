% program for studying SAF-FCS in a spherical confinement

radiusv = [0.1 0.2 0.5 1 2];
rad = (0.5:300)'/300;

maxm = 20; %100;
[weight,x] = gausslegendrecof(2*maxm);
bev = acos(x);
%clear tst; for j=0:maxL for k=0:j tst(j+1,k+1)=sum(weight'*(SphericalHarmonic(j,k,bev,alv).*conj(SphericalHarmonic(j,k,bev,alv)))); end; end
clear M
for j=1:maxm+1
    tmp = legendre(j-1,x,'norm');
    M(j,:)= weight'.*tmp(1,:);
end

maxzeros = 10; %50;
clear kalpha
for j=1:maxm+1
    tmp = mBesselZero(j-1,1e-20,maxzeros+1,3);
    if j>1 && tmp(1)<1e-5
        kalpha(j,:) = tmp(2:end);
    else
        kalpha(j,:) = tmp(1:end-1);
    end
end

clear nrm
for j=1:maxm+1
    for k=1:maxzeros
        if j==1 && k==1
            nrm(j,k) = sqrt(2/3);
        else
            nrm(j,k) = sqrt(pi*besselj(j-0.5,kalpha(j,k))^2/kalpha(j,k)/4*(1-j*(j-1)/kalpha(j,k)^2));
        end
    end
end

clear Jv
dr = mean(diff(rad));
for j=1:maxm+1
    for k=1:maxzeros
        if j==1 && k==1
            Jv(:,k,j) = dr*rad.^2/nrm(j,k);
        else
            Jv(:,k,j) = dr*rad.^2.*sqrt(pi/2./(kalpha(j,k).*rad)).*besselj(j-0.5,kalpha(j,k)*rad)/nrm(j,k);
        end
    end
end

Nsub = 3;
lam = 2^(1/Nsub);
timewin = 1e7;
Ncasc = ceil(log2(timewin));
autotime = [0; lam.^(0:Ncasc*Nsub-1)'];
dif = 5e-7; % diffusion constant in mum^2/time unit
tirf = 0.1;

auto = zeros(length(autotime),length(radiusv));
intensity = zeros(1,length(radiusv));
for jr=1:length(radiusv)
    radius = radiusv(jr); % cell radius in mum
    zv = radius*+rad*x';
    tmp = exp(-zv/tirf);
    
    coef = zeros(maxm+1,maxzeros);
    for j=1:maxm+1
        coef(j,:) = abs(Jv(:,:,j)'*tmp*M(j,:)')';
    end
    
    intensity(jr) = coef(1,1);
    for j=1:maxm+1
        for k=1:maxzeros
            auto(:,jr) = auto(:,jr) + coef(j,k)^2*exp(-dif*(kalpha(j,k)/radius).^2*autotime);
        end
    end
end
semilogx(autotime,auto./(ones(size(auto,1),1)*auto(1,:)))
xlabel('time (a.u.)')
ylabel('autocorrelation');
s={};
for j=1:length(radiusv)
    s{j} = ['\ita\rm = ' num2str(radiusv(j)) ' \mum'];
end
legend(s,3)

figure
plot(radiusv,auto(end,:)./auto(1,:).*radiusv.^3)
grid
xlabel('vesicle radius (µm)');
ylabel('effective detection volume (µm^3)')

figure
semilogy(radiusv,1./(auto(end,:)./auto(1,:).*radiusv.^3*1e-15)/6.022e23)
grid
xlabel('vesicle radius (µm)');
ylabel('single molecule concentration (M)')

figure
semilogy(radiusv,4*pi/3./(auto(end,:)./auto(1,:)))
grid
xlabel('vesicle radius (µm)');
ylabel('molecules per vesicle @ SM concentration')
