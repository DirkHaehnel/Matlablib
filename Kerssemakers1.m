% program for checking DIAM for height measurements 

close all 
clear all

rhov = [0 5];
zv = 0; %(0:5:300)/1e3;
NA = 1.45;
n1 = 1.52;
n = 1.33;
n2 = n;
d1 = [];
d = 0;
d2 = [];
lamem = 0.63;
mag = 60;
focpos = 0;

for j=1:length(zv)
    z = zv(j);
    [tmp, tmp, tmp, rho, phi, fxx0a, fxx2a, fxza, byx0a, byx2a, byza] = SEPDipole(rhov, z, NA, n1, n, n2, d1, d, d2, lamem, mag, focpos);
    % [tmp, tmp, tmp, rho, phi, fxx0b, fxx2b, fxzb, byx0b, byx2b, byzb] = SEPDipole(rhov, z, NA, n1, n, n2, d1, d, d2, lamem, mag, focpos, [], [], [n pi/2]);
    [tmp, tmp, tmp, rho, phi, fxx0b, fxx2b, fxzb, byx0b, byx2b, byzb] = SEPDipole(rhov, z, n, n1, n, n2, d1, d, d2, lamem, mag, focpos);
    inta(j,:) = real(fxx0a.*conj(byx0a) + fxx2a.*conj(byx2a) + fxza.*conj(byza));
%     fxx0a = fxx0a - fxx0b;
%     fxx2a = fxx2a - fxx2b;
%     fxza = fxza - fxzb;
%     byx0a = byx0a - byx0b;
%     byx2a = byx2a - byx2b;
%     byza = byza - byzb;
%     intd(j,:) = real(fxx0a.*conj(byx0a) + fxx2a.*conj(byx2a) + fxza.*conj(byza));
    intb(j,:) = real(fxx0b.*conj(byx0b) + fxx2b.*conj(byx2b) + fxzb.*conj(byzb));
    plot(rho,inta,rho,intb); drawnow;
end

return 

plot(zv,max(inta')'./max(intb')'-1)
plot(zv,(inta*rho')./(intb*rho')-1)

nn = 50;
[x,y] = meshgrid(-nn:nn,-nn:nn);
x = x*8;
y = y*8;
p = angle((x-100)+i*y);
r = sqrt((x-100).^2+y.^2);
rho = rho(1):diff(rho(1:2)):(max(r(:))+diff(rho(1:2)));
inta = [inta zeros(1,length(rho)-size(inta,2))];
int = interp1(rho,inta(j,:),r,'cubic');