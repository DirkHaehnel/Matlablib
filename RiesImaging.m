% program for checking DIAM for height measurements 

close all 
clear all

%NA = 1.49;
NA = 1.65;
%n1 = 1.5;
n1 = 1.78;
n = 1.33;
n2 = n;
d1 = [];
d = 0;
d2 = [];
lambda = 0.67;
mag = 300;
focpos = 0;
atf = [];
pix = 16;
nn = 20;

if 0 % defoc images as function of z-position
    ring = [];
    focpos = 0.5;
    zv = [0:5:100]/1e3;
    for jz = 1:length(zv)
        z = zv(jz);
        intxs = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, 0]);
        intys = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, pi/2]);
        intzs = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [0, 0]);
        CombineImages(cat(3,intxs,intys,intzs,intxs+intys+intzs),1,4,[],{'\itx','\ity','\itz',['iso: \itz\rm = ' mint2str(1e3*z,3,' ') ' nm']});
        eval(['print -dpng -r200 Ries' mint2str(jz,2)]);
    end
end

if 0 % Jonas' problem
    focpos = 0;
    zv = (0:20:100)/1e3;
    for jz = 1:length(zv)
        z = zv(jz);
        ring = [];
        % full
        intf(:,:,jz) = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, 0]);
        intf(:,:,jz) = intf(:,:,jz) + DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, pi/2]);
        intf(:,:,jz) = intf(:,:,jz) + DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [0, 0]);
        % UAF
        intc(:,:,jz) = DefocImage(nn, pix, z, n, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, 0]);
        intc(:,:,jz) = intc(:,:,jz) + DefocImage(nn, pix, z, n, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, pi/2]);
        intc(:,:,jz) = intc(:,:,jz) + DefocImage(nn, pix, z, n, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [0, 0]);
        % ring1
        ring = '1i*1e3*(rad<n/NA)';
        int1(:,:,jz) = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, 0]);
        int1(:,:,jz) = int1(:,:,jz) + DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, pi/2]);
        int1(:,:,jz) = int1(:,:,jz) + DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [0, 0]);
    end
    mm = max(squeeze(sum(sum(intf,1),2)));
    plot(1e3*zv,squeeze(sum(sum(intf,1),2))/mm,1e3*zv,squeeze(sum(sum(intc,1),2))/mm,1e3*zv,squeeze(sum(sum(int1,1),2))/mm)
    xlabel('distance from surface (nm)')
    ylabel('image intensity (a.u.)');
    legend({'all light','UAF','SAF'})
    figure
    CombineImages(cat(3,intf,intc,int1),3,6,[],{'0 nm','20 nm','40 nm','60 nm','80 nm','100 nm'},{'full','UAF','SAF'});
    figure
    CombineImages(cat(3,intf,intc,int1),3,6,'scale',{'0 nm','20 nm','40 nm','60 nm','80 nm','100 nm'},{'full','UAF','SAF'});
end

if 0 % dependence on ring aperture
    focpos = 0;
    zv = (0:20:100)/1e3;
    opening = 0.5:0.1:1;  
    for jz = 1:length(zv)
        z = zv(jz);
        ring = [];
        % full
        int = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, 0]);
        int = int + DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, pi/2]);
        int = int + DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [0, 0]);
        amp_all(jz) = sum(sum(int));
        % UAF
        int = DefocImage(nn, pix, z, n, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, 0]);
        int = int + DefocImage(nn, pix, z, n, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, pi/2]);
        int = int + DefocImage(nn, pix, z, n, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [0, 0]);
        amp_uaf(jz) = sum(sum(int)); 
        for jo = 1:length(opening)
            % ring
            ring = ['1i*1e3*(rad<' num2str(opening(jo)) '*n/NA)'];
            int = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, 0]);
            int = int + DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, pi/2]);
            int = int + DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [0, 0]);
            amp_ring(jz,jo) = sum(sum(int));
        end
    end
    mm = max(amp_all);
    plot(1e3*zv,amp_all./amp_uaf,1e3*zv,amp_ring./(amp_uaf'*ones(1,size(amp_ring,2))));
    colorize
    xlabel('distance from surface (nm)')
    ylabel('rel. image intensity (UAF)');
end

if 1 % dependence on ring aperture: images
    focpos = 0;
    z = 0;
    opening = 0.5:0.1:1;
    NAv = [1.49 1.65];
    n1v = [1.5 1.78];
    for jna=1:length(NAv)
        NA = NAv(jna);
        n1 = n1v(jna);
        ring = [];
        % full
        int = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, 0]);
        int = int + DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, pi/2]);
        int_all(:,:,jna) = int + DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [0, 0]);
        % UAF
        int = DefocImage(nn, pix, z, n, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, 0]);
        int = int + DefocImage(nn, pix, z, n, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, pi/2]);
        int_uaf(:,:,jna) = int + DefocImage(nn, pix, z, n, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [0, 0]);
        for jo = 1:length(opening)
            % ring
            ring = ['1i*1e3*(rad<' num2str(opening(jo)) '*n/NA)'];
            int = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, 0]);
            int = int + DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [pi/2, pi/2]);
            int_ring(:,:,jo,jna) = int + DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, [0, 0]);
        end
    end
    CombineImages(cat(3,int_all(:,:,1),int_uaf(:,:,1),int_ring(:,:,:,1),int_all(:,:,2),int_uaf(:,:,2),int_ring(:,:,:,2)),2,8,...
        [],{{'all',''},{'UAF',''},{'0.5',''},{'0.6',''},{'0.7',''},{'0.8',''},{'0.9',''},{'1.0',''}},{'NA = 1.49','NA = 1.65'});
    % print -dpng -r600 JonasNA149VariableAperturePSF
end


