% program PhasePlateExc

zfield = 0;
NA = 1.3;
fd = 1.8e3;
n0 = 1.51;
n = 1.51;
n1 = 1.51;
d0 = [];
d = 0;
d1 = [];
lamex = 0.633;
% lamex = 0.488;
over = 2.5e3;
focpos = 0;
atf = [];
resolution = [25 5];
rhofield = [0-lamex/resolution(1)/2 1.2];

tshv = (2.9e3:3e3)/NA/fd/2; clear tstx
for j=1:length(tshv)
    tsh = tshv(j);
    ring = ['-i*pi.*(rad<' num2str(tsh) ')'];
    [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc([0-lamex/resolution(1)/2 lamex/resolution(1)/2], zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring);
    tstx(j) = abs(fxc(1,1,1));
end    
    
tsh = mean(tshv(tstx==min(tstx)));
ring = ['-i*pi.*(rad<' num2str(tsh) ')'];
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring);
phi = 0:pi/100:2*pi;   
row = ones(size(phi));

fx = fxc(:,1,1)*row; fy = fyc(:,1,1)*row; fz = fzc(:,1,1)*row;
for j=1:size(fxs,3)
    fx = fx + fxc(:,:,j+1)*cos(j*phi) + fxs(:,1,j)*sin(j*phi);
    fy = fy + fyc(:,:,j+1)*cos(j*phi) + fys(:,1,j)*sin(j*phi);
    fz = fz + fzc(:,:,j+1)*cos(j*phi) + fzs(:,1,j)*sin(j*phi);
end

% bev = (0:30:90)/180*pi;
% alv = (90:-30:0)/180*pi;

mm = max(max(abs(fx).^2));

bev = (0:15:90)/180*pi;
alv = (90:-15:0)/180*pi;
clf
for ja=1:length(alv)
    for jb=1:length(bev)
        %subplot(length(alv),length(bev),length(bev)*(ja-1)+jb)    
        subplot('position',[(jb-1)/length(bev) 1-ja/length(alv) 0.9/length(bev) 0.9/length(alv)])    
        ff = (sin(alv(ja))*cos(bev(jb))*abs(fx)).^2+(sin(alv(ja))*sin(bev(jb))*abs(fy)).^2+(cos(alv(ja))*abs(fz)).^2;
        mpcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),ff);
        text(1.0,-0.9,mnum2str(max(max(ff))/mm,1,2),'FontSize',10,'HorizontalAlignment','left')
        axis off
    end
end

return

bev = (80:2:90)/180*pi;
alv = (90:-5:70)/180*pi;
clf
for ja=1:length(alv)
    for jb=1:length(bev)
        %subplot(length(alv),length(bev),length(bev)*(ja-1)+jb)    
        subplot('position',[(jb-1)/length(bev) 1-ja/length(alv) 0.9/length(bev) 0.9/length(alv)])    
        ff = (sin(alv(ja))*cos(bev(jb))*abs(fx)).^2+(sin(alv(ja))*sin(bev(jb))*abs(fy)).^2+(cos(alv(ja))*abs(fz)).^2;
        mpcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),ff);
        text(1.0,-0.9,mnum2str(max(max(ff))/mm,1,2),'FontSize',12,'HorizontalAlignment','left')
        axis off
    end
end

