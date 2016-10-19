if 0 % Rabeau

    lamex = 0.76;
    resolution = 50;
    rhofield = [-lamex/resolution(1)/2 3.];
    molpos = 15;
    pixel = 0.01; %/sqrt(2);
    NA = 1.3;
    fd = 3e3;
    n0 = 1.51;
    n = 2.5;
    n1 = 2.5;
    d0 = [];
    d = 0;
    d1 = [];
    over = inf;
    focposv = molpos*n0/n - (0:47)*0.05;
    atf = [];
    ring = [];
    psi = (180-28)/180*pi;

else % Olaf Schulz

    lamex = 0.64;
    resolution = 50;
    rhofield = [-lamex/resolution(1)/2 3.];
    molpos = 0;
    pixel = 0.01; %/sqrt(2);
    NA = 1.45;
    fd = 3e3;
    n0 = 1.51;
    n = 1.0;
    n1 = 1.0;
    d0 = [];
    d = 0;
    d1 = [];
    over = inf;
    focposv = -(0:47)*0.02;
    atf = [];
    ring = [];
    psi = 0;
    
end

% circular polarized excitation:
% [fxcr, fxsr, fycr, fysr, fzcr, fzsr] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
% fxc = fxc + i*fxcr;
% fxs = fxs + i*fxsr;
% fyc = fyc + i*fycr;
% fys = fys + i*fysr;
% fzc = fzc + i*fzcr;
% fzs = fzs + i*fzsr;

if 1
    theta = (90:-5:0)/180*pi;
    % in-plane rotation of dipole axis
    phi = 0*pi/2; %(0:5:355)/180*pi;
    theta = ones(length(phi),1)*theta;
    theta = theta(:)';
    phi = repmat(phi,[1 length(theta)/length(phi)]);
    nn = 100; %floor(rhofield(2)/pixel/sqrt(2));
    mask = zeros(2*nn+1,2*nn+1,length(focposv),length(theta));
    for jz=1:length(focposv)
        [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, molpos, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposv(jz), atf, resolution);
        
        for cnt=1:length(theta)
            al = theta(cnt);
            be = phi(cnt);
            mask(:,:,jz,cnt) = DefExcImage(al,be,nn,pixel,rho,fxc,fxs,fyc,fys,fzc,fzs,psi);
        end
    end
    
    %     for j=1:8 sx{j} = [int2str((max(focposv)-focposv(j))*1e3) ' nm']; end
    %     for j=1:6 sy{j} = ['+' int2str((focposv(1)-focposv(9))*1e3) ' nm']; end
    for j=1:8 sx{j} = [int2str(focposv(j)*1e3) ' nm']; end
    for j=1:6 sy{j} = ['+' int2str(focposv(9)*1e3) ' nm']; end
    sy{1} = '';
    CombineImages(mask(:,:,:,11),6,8,'scale',sx,sy);
    
else
    mask = zeros(2*nn+1,2*nn+1,size(z,2),6);
    [int, x, y, feldx, feldy, feldz] = DefExcImage(al,be,nn,pixel,rho,fxc,fxs,fyc,fys,fzc,fzs);
    mask(:,:,:,1) = abs(feldx).^2;
    mask(:,:,:,2) = abs(feldy).^2;
    mask(:,:,:,3) = abs(feldz).^2;
    mask(:,:,:,4) = real(feldx.*conj(feldy) + feldy.*conj(feldx));
    mask(:,:,:,5) = real(feldy.*conj(feldz) + feldz.*conj(feldy));
    mask(:,:,:,6) = real(feldz.*conj(feldx) + feldx.*conj(feldz));
end

if 0
    for k = 1:size(z,2)
        col = ceil(sqrt(size(mask,4)));
        wdth = size(mask,2);
        hght = size(mask,1);
        im = zeros(ceil(size(mask,4)/col)*hght,col*wdth);
        for j=1:size(mask,4)
            im(fix((j-1)/col)*hght+1:(fix((j-1)/col)+1)*hght,mod(j-1,col)*wdth+1:(mod(j-1,col)+1)*wdth) = mask(:,:,k,j)/max(max(mask(:,:,k,j)));
        end
        mim(im)
        eval(['print -dpng -r300 rabeau' mint2str(k,2)]);
    end
end

return

k=28;
tmp=reshape(squeeze(mask(:,:,k,:)),201^2,6);
c=tmp\tst(:);
bla=0*tst; for j=1:6 bla=bla+c(j)*mask(:,:,k,j); end
mim(cat(3,tst,bla))

return

%search patterns
load c:\Joerg\Doc\Microscopy\Rabeau\scans.mat scan_z_1320 scan_z_720 scan_z_n500
im = scan_z_1320;
clear err bim cim sim xc yc bc cc sc len imm
if 1
    k=28;
    [err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm] = FindPattern(im,squeeze(mask(:,:,k,:)),[],[],[],1,[],0);
else
    for k = 1:size(mask,3)
        [err, bim, cim, sim, xc{k}, yc{k}, bc, cc, sc{k}, len, imm{k}] = FindPattern(im,squeeze(mask(:,:,k,:)),[],[],[],1,[],0);
    end
end

return

% check results
for j=1:length(xc) if ~isempty(xc{j}) mim(cat(3,im,imm{j})); title(j); pause; end; end

return

% plot best result
j = 426; k = 28; im = scan_z_1320;
mim(cat(3,im(50:end,20:end-30),mask(:,:,k,j)));
title({['\phi = ' num2str(phi(j)/pi*180) '°, \theta = ' num2str(theta(j)/pi*180) '°, \Delta\itz\rm = ' num2str(z(1,k)) ' \mum']})


j = 426; k = 16; im = scan_z_720;
mim(cat(3,im(50:end,20:end-30),mask(:,:,k,j)));
title({['\phi = ' num2str(phi(j)/pi*180) '°, \theta = ' num2str(theta(j)/pi*180) '°, \Delta\itz\rm = ' num2str(z(1,k)) ' \mum']})


j = 426; k = 11; im = scan_z_n500;
mim(cat(3,im(50:end,20:end-30),mask(:,:,k,j)));
title({['\phi = ' num2str(phi(j)/pi*180) '°, \theta = ' num2str(theta(j)/pi*180) '°, \Delta\itz\rm = ' num2str(z(1,k)) ' \mum']})
