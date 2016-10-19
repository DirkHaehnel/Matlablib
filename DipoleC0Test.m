% Testing the calculation of the angular distribution of radiation of a
% dipole in clyindrical vector functions

clear all
close all

n = 1.33;
rho = 1500/570*2*pi;
theta = (pi/500:pi/250:pi)';
st = sin(theta);
ct2 = cos(theta).^2;
col = ones(size(theta));
nmax = 30;
w = n(end)*cos(theta);
q = sqrt(n^2 - w.^2);

if 0
    for m=-nmax:nmax
        c0(:,nmax+m+1,1) = CylM(m,q,rho,'j')'./q'.^2;
        c0(:,nmax+m+1,2) = CylN(m,q,w,rho,'j')'./q'.^2/n^2;
    end
end

if 1
    for m=-nmax:nmax
        c0(:,nmax+m+1,1) = CylM(m,q,rho,'j','f')'./q'.^2;
        c0(:,nmax+m+1,2) = CylN(m,q,w,rho,'j','f')'./q'.^2/n^2;
    end
end

if 0
    for m=-nmax:nmax
        c0(:,nmax+m+1,1) = 0*q;
        c0(:,nmax+m+1,2) = CylN(m,q,w,rho,'j','z')'./q'.^2/n^2;
    end
end

phi = pi/500:pi/250:2*pi;
% er = zeros(length(theta),length(phi),2);
% for j=1:2*nmax+1
%     er(:,:,1) = er(:,:,1) + (st.*c0(:,j,1))*exp(i*(j-nmax-1)*(phi-pi/2));
%     er(:,:,2) = er(:,:,2) + n*(st.*c0(:,j,2))*exp(i*(j-nmax-1)*(phi-pi/2));    
% end
% er = n/pi*er;
% adr = pi/2*n*squeeze((abs(er(:,:,1)).^2 + abs(er(:,:,2)).^2));
% adr = Cyl2ADR(theta,phi,er,mv,n)
adr = Cyl2ADR(theta,phi,c0,-nmax:nmax,n);

surf((sin(theta)*cos(phi)).*adr, (sin(theta)*sin(phi)).*adr, (cos(theta)*ones(1,length(phi))).*adr,adr);
shading interp
camlight
axis image
cameratoolbar

tst=sum(sin(theta)'*adr)*mean(diff(theta))*mean(diff(phi))

