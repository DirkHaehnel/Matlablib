% Program PDS for ab initio modeling of Photon Distance Statistics data
close all

focpos = 0; % Fokus-Position
over = 0.99; % aperture filling
n0 = 1.33; n1 = 1.33;
NA = 1.2;
lam = 0.635;
%[exc, rho, z] = hecht([0 4*lam], [0 4*lam], NA, n0, n1, lam, over, focpos); 

dif = 1e-4; % diffusion constant in mum^2/time unit
av = [50 100 200 300]/66.667; % aperture diameters
f0 = exc;
f0 = repmat(f0,[1 1 length(av)]);
rhom = repmat(rho,[1 1 length(av)]);
for j=1:length(av)
   f0(:,:,j) = d3cef(rho/lam*2*pi,z/lam*2*pi,av(j)/lam*2*pi,NA).*f0(:,:,j);
   pcolor(rho,z,squeeze(f0(:,:,j))); axis image; shading interp; drawnow; 
end

Nsub = 3;
lam = 2^(1/Nsub);
timewin = 1e7;
Ncasc = ceil(log2(timewin));
autotime = lam.^(0:Ncasc*Nsub-1)';
close
Ncoef = 10;
pdscoef = zeros(length(autotime),length(av),Ncoef);
dz = z(1,2)-z(1,1); drho = rho(2,1)-rho(1,1);
for j = 1:length(autotime)
   
   tmp = (erf((z-z'+dz/2)/sqrt(4*dif*autotime(j))) - erf((z-z'-dz/2)/sqrt(4*dif*autotime(j))))/2;
   tmp = tmp + (erf((z+z'+dz/2)/sqrt(4*dif*autotime(j))) - erf((z+z'-dz/2)/sqrt(4*dif*autotime(j))))/2;
   for k1 = 1:size(f0,1)
      for k2 = 1:size(f0,3)
         f(k1,:,k2) = f0(k1,:,k2)*tmp;
      end
   end
   
   rr = rho(:,1)*rho(:,1)'/2/dif/autotime(j);
   tmp = (rho.^2+rho'.^2)/4/dif/autotime(j);
   tmp1 = tmp - rr;
   tmp(rr<=700) = besselj(0, i*rr(rr<=700)).*exp(-tmp(rr<=700));
   tmp(rr>700) = exp(-tmp1(rr>700))./sqrt(rr(rr>700)*2*pi);
   tmp = (ones(size(rho(:,1)))*rho(:,1)').*tmp/2/dif/autotime(j)*drho;
   ss = sum(tmp')';
   tmp(ss>1,:) = diag(1./ss(ss>1))*tmp(ss>1,:); 
   for k1 = 1:size(f0,2)
      for k2 = 1:size(f0,3)
         f(:,k1,k2) = tmp*f(:,k1,k2);
      end
   end
   for k = 1:Ncoef
      pdscoef(j,:,k) = (-1)^k*squeeze(sum(sum(rhom.*f.^k)))';
   end
   j 
   squeeze(pdscoef(j,:,:))
end
pdscoef = pdscoef*2*pi*drho*dz;

%tmp=ones(size(autotime))*(1./cumprod(1:size(pdscoef,3)));
%tmp1=ones(size(autotime))*(1./cumprod([1 1:size(pdscoef,3)-1]));
%loglog(autotime,-fac*sum((tmp1.*squeeze(pdscoef(:,1,:)))').*exp(fac*sum((tmp.*squeeze(pdscoef(:,1,:)))')))

%save d:/jender/temp/foci100 exc pdscoef dif rho z lam n0 n1 NA drho av rhov zv

