% program ExcitonMicroscopy

close all
bld = 1;

k10 = [1/20 1/0.7 1/0.5];
extinction = 2.5e5*[1 1 1];
Tpulse = 0.05;
Trep = Tpulse; 
%Trep = 12.5;
Imax = 10;

rhofield = [0 1];
zfield = [8 12];
NA = 1.2;
fd = 3e3;
n0 = 1.33;
n = n0;
n1 = n0;
d0 = []; 
d = 0;
d1 = [];
lamex = 0.47;
over = 5e3;
focpos = 10;
exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, [], 50);

rho = exc.rho; 
z = exc.z-focpos;
[feld, phi, rr, zz] = FocusImage3D(rho,z,cat(4,cat(3,exc.fxc,exc.fxs),cat(3,exc.fyc,exc.fys),cat(3,exc.fzc,exc.fzs)));
feld = sum(abs(feld).^2,4);
feld = feld/max(feld(:));

% close; FocusImage3D(exc.rho,exc.z,cat(4,cat(3,exc.fxc,exc.fxs),cat(3,exc.fyc,exc.fys),cat(3,exc.fzc,exc.fzs)))

exc.fxc(:,:,1) = 1;
exc.fxc(:,:,2:end) = 0;
exc.fyc(:) = 0;
exc.fzc(:) = 0;
exc.fxs(:) = 0;
exc.fys(:) = 0;
exc.fzs(:) = 0;

lamem = 0.6;
mag = 60;
av = 25;
mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av);
cef = FocusImage3D(rho,z,mdf.volx+mdf.voly);
cef = cef/max(cef(:));

% close; FocusImage3D(exc.rho,exc.z,mdf.volx+mdf.voly)

v = Imax*10.^(-10:0.01:0);
[Ss If] = ExcitonExcitation(v,k10,extinction,Tpulse,Trep);

if bld
    plot(v,If)
    xlabel('excitation intensity (\muW/\mum^2)');
    ylabel('photon emissions (1/s)')
    legend(2,{'exciton','bi-exciton','tri-exciton'})
end

if bld
    figure
    ind = sum(v<1e-2);
    loglog(v(ind:end),If(:,ind:end))
    hold on
    loglog(v(ind:end),If(1,ind)/v(ind)*v(ind:end),':r');
    loglog(v(ind:end),If(2,ind)/v(ind)^2*v(ind:end).^2,':b');
    loglog(v(ind:end),If(3,ind)/v(ind)^3*v(ind:end).^3,':g');
    hold off
    ax = axis;
    axis([1e-2 1e1 10^floor(log10(min(If(:,ind)))) ax(4)])
    xlabel('excitation intensity (\muW/\mum^2)');
    ylabel('photon emissions (1/s)')
    legend(4,{'exciton','bi-exciton','tri-exciton'})
end

Iexc = 0.0681; % corresponds to ca. 1 muW/mum^2
f1 = interp1(v,If(1,:),Iexc*feld,'cubic');
f1 = f1/max(f1(:));
f2 = interp1(v,If(2,:),Iexc*feld,'cubic');
f2 = f2/max(f2(:));
f3 = interp1(v,If(3,:),Iexc*feld,'cubic');
f3 = f3/max(f3(:));

if bld
   figure
   FocusImage3D(rr,zz,feld.*cef,phi); 
   ax = axis;
   figure
   FocusImage3D(rr,zz,f3.*cef,phi); 
   axis(ax) 
end

x = [-squeeze(rr(end:-1:1,round(size(rr,2)/2+1),1)); squeeze(rr(:,round(size(rr,2)/2+1),1))];
tmp0 = squeeze(feld(:,round(size(rr,2)/2+1),1).*cef(:,round(size(rr,2)/2+1),1));
tmp0 = [tmp0(end:-1:1); tmp0];
tmp1 = squeeze(f1(:,round(size(rr,2)/2+1),1).*cef(:,round(size(rr,2)/2+1),1)); 
tmp1 = [tmp1(end:-1:1); tmp1];
tmp2 = squeeze(f2(:,round(size(rr,2)/2+1),1).*cef(:,round(size(rr,2)/2+1),1)); 
tmp2 = [tmp2(end:-1:1); tmp2];
tmp3 = squeeze(f3(:,round(size(rr,2)/2+1),1).*cef(:,round(size(rr,2)/2+1),1)); 
tmp3 = [tmp3(end:-1:1); tmp3];
wx0 = sqrt(sum(x.^2.*tmp0)/sum(tmp0)); 

if bld
    figure
    plot(x,tmp0,x,tmp2,x,tmp3)
    axis([-0.5 0.5 0 1])
    xlabel('lateral position (\mum)');
    ylabel('rel. fluorescence intensity')
end

y = squeeze(zz(1,:,1));
tmp0 = squeeze(feld(1,:,1).*cef(1,:,1));
tmp1 = squeeze(f1(1,:,1).*cef(1,:,1)); 
tmp2 = squeeze(f2(1,:,1).*cef(1,:,1)); 
tmp3 = squeeze(f3(1,:,1).*cef(1,:,1)); 
wz0 = sqrt(sum(y.^2.*tmp0)/sum(tmp0)); 

if bld
    figure
    plot(y,tmp0,y,tmp2,y,tmp3)
    xlabel('position along optical axis (\mum)');
    ylabel('rel. fluorescence intensity')
end

for j=1:100
    Iexc = j/100*Imax;
    f = interp1(v,If(1,:),Iexc*feld,'cubic');
    f = f/max(f(:));
    tmp = squeeze(f(:,round(size(rr,2)/2+1),1).*cef(:,round(size(rr,2)/2+1),1));
    tmp = [tmp(end:-1:1); tmp];
    wx1(j) = sqrt(sum(x.^2.*tmp)/sum(tmp));
    tmp = squeeze(f(1,:,1).*cef(1,:,1));
    wz1(j) = sqrt(sum(y.^2.*tmp)/sum(tmp));
    f = interp1(v,If(2,:),Iexc*feld,'cubic');
    f = f/max(f(:));
    tmp = squeeze(f(:,round(size(rr,2)/2+1),1).*cef(:,round(size(rr,2)/2+1),1));
    tmp = [tmp(end:-1:1); tmp];
    wx2(j) = sqrt(sum(x.^2.*tmp)/sum(tmp));
    tmp = squeeze(f(1,:,1).*cef(1,:,1));
    wz2(j) = sqrt(sum(y.^2.*tmp)/sum(tmp));
    f = interp1(v,If(3,:),Iexc*feld,'cubic');
    f = f/max(f(:));
    tmp = squeeze(f(:,round(size(rr,2)/2+1),1).*cef(:,round(size(rr,2)/2+1),1));
    tmp = [tmp(end:-1:1); tmp];
    wx3(j) = sqrt(sum(x.^2.*tmp)/sum(tmp));
    tmp = squeeze(f(1,:,1).*cef(1,:,1));
    wz3(j) = sqrt(sum(y.^2.*tmp)/sum(tmp));
end

if bld
    figure
    plot((1:100)/100*Imax*pi*wx0^2*2, 2*wx1,(1:100)/100*Imax*pi*wx0^2*2, 2*wx2,(1:100)/100*Imax*pi*wx0^2*2, 2*wx3)
    xlabel('excitation power (\muW)');
    ylabel('lateral resolution (\mum)')
    grid
    figure
    plot((1:100)/100*Imax*pi*wx0^2*2, 2*wz1,(1:100)/100*Imax*pi*wx0^2*2, 2*wz2,(1:100)/100*Imax*pi*wx0^2*2, 2*wz3)
    xlabel('excitation power (\muW)');
    ylabel('axial resolution (\mum)')
    grid
end

return

[y z] = ExcitonExcitation(x,k10,extinction,Tpulse,Trep)

close; subplot(121); FocusImage3D(rr,zz,feld.*cef,phi,[],[],1); subplot(122); FocusImage3D(rr,zz,tst.*cef,phi,[],[],1); axis(ax)