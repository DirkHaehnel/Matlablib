function Field2Image(rho,z,fc,fs,tsh,varargin);

if nargin<5 | isempty(tsh)
    tsh = 1/exp(2);
end

if isempty(fc) 
    fc = zeros(size(fs,1),size(fs,2),1);
elseif isempty(fs) 
    fs = zeros(size(fc,1),size(fc,2),1);
end
    
tst = sum(abs(fc),3).^2+sum(abs(fs),3).^2;
tst = tst/max(tst(:));

dz = z(1,2)-z(1,1);
rhomax = 2*max(rho(tst>=tsh));
zmax = max(z(tst>=tsh))+10*dz;
zmin = min(z(tst>=tsh))-10*dz;

indz = z(1,:)>=zmin & z(1,:)<=zmax; indr = rho(:,1)<=rhomax;
rho = rho(indr,indz);
z = z(indr,indz);
fc = fc(indr,indz,:);
fs = fs(indr,indz,:);

rr = rho(rho(:,1)<rho(end,1)/sqrt(2),1); 
rr = [-flipud(rr); rr];
[xx,yy,zz] = meshgrid(rr,rr,z(1,:)');
rr = sqrt(xx.^2+yy.^2); phi = angle(xx+1i*yy);

feld = interp2(z,rho,fc(:,:,1),zz,rr,'cubic');
for j=2:size(fc,3) 
    feld = feld + interp2(z,rho,fc(:,:,j),zz,rr,'cubic').*cos((j-1)*phi); 
end
for j=1:size(fs,3) 
    feld = feld + interp2(z,rho,fs(:,:,j),zz,rr,'cubic').*sin(j*phi); 
end
feld = abs(feld).^2;
%feld = real(feld);

[tx, ty, tz] = surfcv(xx,yy,zz,feld,tsh*max(feld(:)));
fill3(tx, ty, tz, sqrt(tx.^2+ty.^2))
set(gcf,'renderer','zbuffer');    
axis image
box
shading interp;
camlight
lighting gouraud
times16
axis vis3d
if nargin>4 set(gca,varargin{:}); end
