function [feld, xx, yy] = FocusSection(rho,vol)

tmp = floor(sqrt(max(rho(:,1))^2/2)/mean(diff(rho(:,1))));
[xx, yy] = meshgrid((-tmp:tmp)*mean(diff(rho(:,1))),(-tmp:tmp)*mean(diff(rho(:,1))));

rad = sqrt(xx.^2+yy.^2);
phi = unwrap(angle(xx + 1i*yy));

vi = interp1(rho(:,1), vol, rad, 'cubic', 'extrap');

feld = vi(:,:,:,1);
for jz=1:size(vol,2)
    for j=1:(size(vol,3)-1)/2
        feld(:,:,jz) = feld(:,:,jz) + vi(:,:,jz,j+1).*cos(j*phi) + vi(:,:,:,(end+1)/2+j).*sin(j*phi);
    end
end


return
clear m1 m2 wx1 wx2 wy1 wy2
for j=1:size(mdf1.volx,2) 
    [feld,xx,yy]=FocusSection(mdf1.rho,mdf1.volx(:,j,:)); 
    [m1(j), ~, wx1(j), wy1(j)]=Gauss2D(xx(1,:),yy(:,1),feld);
    [feld,xx,yy]=FocusSection(mdf1.rho,mdf2.voly(:,j,:)); 
    [m2(j), ~, wx2(j), wy2(j)]=Gauss2D(xx(1,:),yy(:,1),feld);
end
plot(mdf1.z(1,:),m1-m2,mdf1.z(1,:),wx1+wy1,mdf1.z(1,:),wx2+wy2)


