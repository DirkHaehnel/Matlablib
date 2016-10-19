% BenScanRead

if 1 % new stuff
    %load c:\Joerg\Doc\Bowen\SurfaceConcentration\joerg_data.mat

    %jd = 30;
    %im = squeeze(mean(double(Data{jd}.Image_Stack(5:end,:,:))));

    im = tiffread('c:\Joerg\Doc\Bowen\SurfaceConcentration\single molecule imaging june 7 2007 large range1.tif');
    im = double(im.data);
    
    [x,y] = meshgrid(-size(im,2)/2+0.5:size(im,2)/2,-size(im,1)/2+0.5:size(im,1)/2);
    smooth = 10;
    im = im - abs(ifft2(ifftshift(fftshift(fft2(im)).*exp(-(x.^2+y.^2)/smooth))));
    im = im - min(im(:));

    fcs2 = mConv2(im,im)/numel(im);
    mm = 20;
    fcs2 = fcs2(round(end/2-mm:end/2+mm),round(end/2-mm:end/2+mm));
    surf(-mm:mm,-mm:mm,fcs2);
    xlabel('\Delta\itx\rm (pix)'); ylabel('\Delta\ity\rm (pix)'); zlabel('\langle\itn\rm(0,0)\cdot\itn\rm(\Delta\itx\rm,\Delta\ity\rm)\rangle')

    ord = 1;
    tmp = fcs2;
    tmp(tmp==max(tmp(:))) = inf;
    close;
    p2 = Simplex('GaussBen',[0 0 2 2 5 5],[-inf -inf 0 0 0 0],[],[],[],tmp,ord);
    p2 = Simplex('GaussBen',p2,[-inf -inf 0 0 0 0],[],[],[],tmp,ord);
    [err, c2, zz2] = GaussBen(p2, tmp, ord);
    xlabel('\Delta\itx\rm (pix)'); ylabel('\Delta\ity\rm (pix)'); zlabel('rel. deviation')

    d = round(3*sqrt(prod(p2(3:4))));
    [x,y] = meshgrid(-d:d,-d:d);
    mask = exp(-(x.^2+y.^2)/prod(p2(3:4)));
    mask((x.^2+y.^2)>d^2) = 0;
    
    for jd = 1:length(Data)
        im = squeeze(mean(double(Data{jd}.Image_Stack(5:end,:,:))));

        [x,y] = meshgrid(-size(im,2)/2+0.5:size(im,2)/2,-size(im,1)/2+0.5:size(im,1)/2);
        smooth = 10;
        im = im - abs(ifft2(ifftshift(fftshift(fft2(im)).*exp(-(x.^2+y.^2)/smooth))));
        im = im - min(im(:));

        [tmp, bim, cim, tmp, resx, resy] = FindPattern(im, mask, ones(size(mask)), 1, 0.5);

        conc(jd) = length(resx)/numel(im);
    end
end

% old variant
if 0
    %name = 'D:\Doc\Bowen\SurfaceConcentration\30sec 2mW 3000 gain dryish 5x in focus.tif';
    %name = 'D:\Doc\Bowen\SurfaceConcentration\e1a.tif';
    %name = 'D:\Doc\Bowen\SurfaceConcentration\e9a.tif';
    name = 'D:\Doc\Bowen\SurfaceConcentration\single molecule imaging june 8th 2007.tif';
    %name = 'D:\Doc\Bowen\SurfaceConcentration\single molecule imaging june 7 2007 large range1.tif';

    im = double(imread(name));
    [x,y] = meshgrid(-size(im,2)/2+0.5:size(im,2)/2,-size(im,1)/2+0.5:size(im,1)/2);
    smooth = 10;
    im = im - abs(ifft2(ifftshift(fftshift(fft2(im)).*exp(-(x.^2+y.^2)/smooth))));
    im = im - min(im(:));

    fcs2 = mconv2(im,im)/numel(im);
    fcs3 = mconv2(im.^2,im)/numel(im);

    mm = 20;
    fcs2 = fcs2(round(end/2-mm:end/2+mm),round(end/2-mm:end/2+mm));
    surf(-mm:mm,-mm:mm,fcs2);
    xlabel('\Delta\itx\rm (pix)'); ylabel('\Delta\ity\rm (pix)'); zlabel('\langle\itn\rm(0,0)\cdot\itn\rm(\Delta\itx\rm,\Delta\ity\rm)\rangle')

    fcs3 = fcs3(round(end/2-mm:end/2+mm),round(end/2-mm:end/2+mm));
    surf(-mm:mm,-mm:mm,fcs3);
    xlabel('\Delta\itx\rm (pix)'); ylabel('\Delta\ity\rm (pix)'); zlabel('\langle\itn\rm(0,0)^2\cdot\itn\rm(\Delta\itx\rm,\Delta\ity\rm)\rangle')

    ord = 1;
    tmp = fcs2;
    tmp(tmp==max(tmp(:))) = inf;
    close;
    p2 = Simplex('GaussBen',[0 0 2 2 5 5],[-inf -inf 0 0 0 0],[],[],[],tmp,ord);
    p2 = Simplex('GaussBen',p2,[-inf -inf 0 0 0 0],[],[],[],tmp,ord);
    [err, c2, zz2] = GaussBen(p2, tmp, ord);
    xlabel('\Delta\itx\rm (pix)'); ylabel('\Delta\ity\rm (pix)'); zlabel('rel. deviation')
    s11 = 2*pi*(prod(reshape(p2(3:end),2,length(p2)/2-1))*c2(end-length(p2)/2+2:end));
    z11 = sum(c2(end-length(p2)/2+2:end));

    tmp = fcs3;
    tmp(tmp==max(tmp(:))) = inf;
    close;
    p3 = Simplex('GaussBen',[0 0 1 1 3 3],[-inf -inf 0 0 0 0],[],[],[],tmp,ord);
    p3 = Simplex('GaussBen',p3,[-inf -inf 0 0 0 0],[],[],[],tmp,ord);
    [err, c3] = GaussBen(p3, tmp, ord);
    xlabel('\Delta\itx\rm (pix)'); ylabel('\Delta\ity\rm (pix)'); zlabel('rel. deviation')
    s21 = 2*pi*(prod(reshape(p3(3:end),2,length(p3)/2-1))*c3(end-length(p3)/2+2:end)) - (2*sqrt(c2(1))+1)*s11;

    d = round(3*sqrt(prod(p2(3:4))));
    [x,y] = meshgrid(-d:d,-d:d);
    mask = exp(-(x.^2+y.^2)/prod(p2(3:4)));
    [tmp, tmp, cim, tmp, resx, resy] = FindPattern(im, mask, ones(size(mask)), 1, 0.5);

    clear tmp; for j=1:length(resx) tmp(j) = cim(resy(j),resx(j)); end
    [h,bin] = mhist(tmp,(0.5:1e2)/1e2*max(tmp));

    conc0 = length(resx)/numel(im);

    tmp = 1:length(bin);
    tmp = tmp(h==max(h)); tmp = tmp(1);
    close; p = simplex('ExpFun',1e3,[],[],[],[],bin(tmp(1):end),h(tmp(1):end),1)
    [err, c] = ExpFun(p, bin(tmp(1):end), h(tmp(1):end));
    hh = h;
    hh(1:tmp(1)) = c(1)+c(2)*exp(-(bin(1:tmp(1))-bin(tmp(1)))/p);
    bar(bin,h); hold on; plot(bin,c(1)+exp(-(bin-bin(tmp(1)))/p)*c(2)); ax=axis; plot(bin(tmp(1))*[1 1],ax(3:4),'g--'); hold off
    xlabel('brightness per molecule (counts)'); ylabel('frequency')
    fac = (sum(bin.^2*hh)^3/sum(bin.^3*hh).^2/sum(hh));
    conc1 = z11^2*s11/s21^2/fac;
end

return

mim(im,'v'); hold on; plot(resx,resy,'.c','markersize',4); hold off

close; p=simplex('ExpFun',1e3,[],[],[],[],bin(3:end),h(3:end),1)
[err, c, zz, z] = ExpFun(p, bin(3:end), h(3:end));
plot(bin,h,bin,c(1)+exp(-(bin-bin(3))/p)*c(2));

return

for k1=-mm:mm
    for k2=-mm:mm
        fcs2(mm+1+k1,mm+1+k2) = sum(sum(im(mm+1+k1:end-mm+k1,mm+1+k2:end-mm+k2).*im(mm+1:end-mm,mm+1:end-mm)));
        fcs3(mm+1+k1,mm+1+k2) = sum(sum(im(mm+1+k1:end-mm+k1,mm+1+k2:end-mm+k2).^2.*im(mm+1:end-mm,mm+1:end-mm)));
    end;
end
fcs2 = fcs2/(size(im,1)-2*mm)/(size(im,2)-2*mm);
fcs3 = fcs3/(size(im,1)-2*mm)/(size(im,2)-2*mm);

return

clear conc fac
for j=1:5
    [imm, cim, resx, resy] = Iris(name, j);
    resx = resx{1}; resy = resy{1};
    conc(j) = length(resx)/numel(imm);
    clear tmp; for k=1:length(resx) tmp(k) = cim(resy(k),resx(k)); end
    fac(j) = sum(tmp.^2)^3/sum(tmp.^3).^2/length(tmp);
end
