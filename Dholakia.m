clear all

maxnum_modes = 40;

rhofield = [0 1.5];
zfield = [-6 6];
NA = 1.2;
fd = 3e3;
n0 = 1.33;
n = n0;
n1 = n0;
d0 = [];
d = 0;
d1 = [];
lamex = 0.64;
over = 3.6e3;
focpos = 0;
atf = [];
resolution = 50;
maxm = 2; % must always correspond to the actually existent modes !!

for j=1:maxnum_modes % calculate electric fields
    s = ['polyval([' num2str(LaguerrePoly(j-1)) '],rad.^2*' num2str(2*(fd*NA/over)^2) ')'];
    s = ['cos(psi).*(' s ')'];
    
    [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, s, maxm);
    [fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
    
    dxc(:,:,:,j) = fxc + fxc2;
    dxs(:,:,:,j) = fxs + fxs2;
    dyc(:,:,:,j) = fyc + fyc2;
    dys(:,:,:,j) = fys + fys2;
    dzc(:,:,:,j) = fzc + fzc2;
    dzs(:,:,:,j) = fzs + fzs2;
end
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, '', maxm);

if 0
    x = 0:0.01:1; 
    j = 1;
    plot(x,polyval(LaguerrePoly(j-1),x.^2*2*(fd*NA/over)^2));
    hold on
    for j=2:maxnum_modes
        plot(x,polyval(LaguerrePoly(j-1),x.^2*2*(fd*NA/over)^2));
    end
    hold off
    colorize
end

if 0
    for j=1:maxnum_modes
        subplot(5,4,j)
        FocusImage2D(rho,z,cat(4,cat(3,dxc(:,:,:,j),dxs(:,:,:,j)),cat(3,dyc(:,:,:,j),dys(:,:,:,j)),cat(3,dzc(:,:,:,j),dzs(:,:,:,j))));
    end
end

for j1=1:maxnum_modes % calculate intensity matrix
    for j2=1:maxnum_modes
        frr(:,:,j1,j2) = real(dxc(:,:,2,j1).*conj(dxc(:,:,2,j2)));
        fzz(:,:,j1,j2) = real(dzc(:,:,1,j1).*conj(dzc(:,:,1,j2)));
    end
end

% calculate goal function matrix
mask = (rho/0.5).^2 + (z/0.1).^2 <= 1;
for j1=1:maxnum_modes 
    for j2=1:maxnum_modes
        m0(j1,j2) = sum(sum(mask(:,end/2).*rho(:,1).*frr(:,end/2,j1,j2)));
    end
end

[v0,e0] = eig(m0);

e0 = diag(e0);
tst = e0>1e-6*max(e0);
ind = 1:length(e0);
ind = ind(tst);

for kk=1:length(ind)
    k = ind(kk);
    j = 1;
    sxc(:,:,:,kk) = v0(j,k)*dxc(:,:,:,j);
    sxs(:,:,:,kk) = v0(j,k)*dxs(:,:,:,j);
    syc(:,:,:,kk) = v0(j,k)*dyc(:,:,:,j);
    sys(:,:,:,kk) = v0(j,k)*dys(:,:,:,j);
    szc(:,:,:,kk) = v0(j,k)*dzc(:,:,:,j);
    szs(:,:,:,kk) = v0(j,k)*dzs(:,:,:,j);
    for j=2:maxnum_modes % calculate electric fields
        sxc(:,:,:,kk) = sxc(:,:,:,kk) + v0(j,k)*dxc(:,:,:,j);
        sxs(:,:,:,kk) = sxs(:,:,:,kk) + v0(j,k)*dxs(:,:,:,j);
        syc(:,:,:,kk) = syc(:,:,:,kk) + v0(j,k)*dyc(:,:,:,j);
        sys(:,:,:,kk) = sys(:,:,:,kk) + v0(j,k)*dys(:,:,:,j);
        szc(:,:,:,kk) = szc(:,:,:,kk) + v0(j,k)*dzc(:,:,:,j);
        szs(:,:,:,kk) = szs(:,:,:,kk) + v0(j,k)*dzs(:,:,:,j);
    end
    sxc(:,:,:,kk) = sxc(:,:,:,kk)/sqrt(e0(k));
    sxs(:,:,:,kk) = sxs(:,:,:,kk)/sqrt(e0(k));
    syc(:,:,:,kk) = syc(:,:,:,kk)/sqrt(e0(k));
    sys(:,:,:,kk) = sys(:,:,:,kk)/sqrt(e0(k));
    szc(:,:,:,kk) = szc(:,:,:,kk)/sqrt(e0(k));
    szs(:,:,:,kk) = szs(:,:,:,kk)/sqrt(e0(k));
end

for j1=1:length(ind) % calculate intensity matrix
    for j2=1:length(ind)
        grr(:,:,j1,j2) = real(sxc(:,:,2,j1).*conj(sxc(:,:,2,j2)));
        gzz(:,:,j1,j2) = real(szc(:,:,1,j1).*conj(szc(:,:,1,j2)));
    end
end

clear m2
for j1=1:length(ind)
    for j2=1:length(ind)
        m2(j1,j2) = sum(sum(mask(:,end/2).*rho(:,1).^3.*(grr(:,end/2,j1,j2)+gzz(:,end/2,j1,j2))));
        %m2(j1,j2) = sum(sum(mask.*(rho.^3.*(grr(:,:,j1,j2)+gzz(:,:,j1,j2)))));
        %m2(j1,j2) = sum(sum(mask.*((rho.^2+(z-focpos).^2)).*rho.*(grr(:,:,j1,j2)+gzz(:,:,j1,j2))));
    end
end

[v2,e2] = eig(m2);

j = 1;
txc = v2(j,1)*sxc(:,:,:,j);
txs = v2(j,1)*sxs(:,:,:,j);
tyc = v2(j,1)*syc(:,:,:,j);
tys = v2(j,1)*sys(:,:,:,j);
tzc = v2(j,1)*szc(:,:,:,j);
tzs = v2(j,1)*szs(:,:,:,j);
for j=2:length(ind) % calculate electric fields
    txc = txc + v2(j,1)*sxc(:,:,:,j);
    txs = txs + v2(j,1)*sxs(:,:,:,j);
    tyc = tyc + v2(j,1)*syc(:,:,:,j);
    tys = tys + v2(j,1)*sys(:,:,:,j);
    tzc = tzc + v2(j,1)*szc(:,:,:,j);
    tzs = tzs + v2(j,1)*szs(:,:,:,j);
end

[sum(sum(mask.*(abs(txc(:,:,2)).^2+abs(tzc(:,:,1)).^2))) sum(sum(mask.*(rho.^2+(z-focpos).^2).*(abs(txc(:,:,2)).^2+abs(tzc(:,:,1)).^2)))/sum(sum(mask.*(abs(txc(:,:,2)).^2+abs(tzc(:,:,1)).^2)))]
subplot(121)
FocusImage2D(rho,z,cat(4,cat(3,txc,txs),cat(3,tyc,tys),cat(3,tzc,tzs)));
colorbar
subplot(122)
FocusImage2D(rho,z,cat(4,cat(3,fxc,fxs),cat(3,fyc,fys),cat(3,fzc,fzs)));
colorbar

p0 = FocusImage2D(rho,z,cat(4,cat(3,txc,txs),cat(3,tyc,tys),cat(3,tzc,tzs)));
p2 = FocusImage2D(rho,z,cat(4,cat(3,fxc,fxs),cat(3,fyc,fys),cat(3,fzc,fzs)));
figure
plot(rho(:,1),p0(:,end/2)/p0(1,end/2),rho(:,1),p2(:,end/2)/p2(1,end/2))

return


close all
FocusImage3D(rho,z,cat(4,cat(3,fxc,fxs),cat(3,fyc,fys),cat(3,fzc,fzs)));
axis([-rhofield(2) rhofield(2) -rhofield(2) rhofield(2) zfield])

return

FocusImage3D(exc.rho,exc.z,cat(4,0*cat(3,exc.fxc,exc.fxs),0*cat(3,exc.fyc,exc.fys),cat(3,exc.fzc,exc.fzs)));

