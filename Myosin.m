function [res, z, zz, a, b, model] = Myosin(file,pic,roi,model,mag,NA)

close all

NA = 1.24; % numerical aperture
n0 = 1.51; % immersion medium
n = 1.0; % sample medium
n1 = n; 
d0 = []; 
d = 0; d1 = [];
lamex = 0.575; % emission wavelength
mag = 110; % magnification
focus = 0.9; % defcousing in mum
pixel = 6.45; % CCD pixel size in mum

be_res = 10; % minimum resolution of in-plane angle
al_res = 10; % minimum resolution of out-of-plane angle

if ~iscell(file)
    filename{1} = file;
else
    filename = file;
end

if nargin<2 | isempty(pic)
    pic = true;
end

if nargin<4 | isempty(model)
    nn = 15; % half size of mask size for pattern matching
    bck = Disk(nn);
    cnt = 1;
    for k=90:-al_res:0
        al = k/180*pi;
        if k==90
            jj = round(180/be_res);
            dbe = pi/jj;
        elseif k==0
            jj = 1;
            dbe = 0;
        else
            jj = round(sin(al)*360/be_res);
            dbe = 2*pi/jj;
        end
        for j=1:jj
            theta(cnt) = al;
            phi(cnt) = dbe*(j-1);
            cnt = cnt+1;
        end
    end

    [intx inty intz, rho, tmp, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 ceil(nn*sqrt(2)*pixel/mag)], 0, NA, n0, n, n1, d0, d, d1, lamex, mag, focus);
    for cnt=1:length(theta)
        al = theta(cnt);
        be = -phi(cnt);        
        mask(:,:,cnt) = SEPImage(al,be,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
        cnt
    end
    if 0
        col = ceil(sqrt(size(mask,3)));
        for j=1:size(mask,3) 
            subplot(col,col,j); 
            mim(mask(:,:,j)); 
            line([0.5 size(mask,2)+0.5 size(mask,2)+0.5 0.5 0.5],[0.5 0.5 size(mask,1)+0.5 size(mask,1)+0.5 0.5],'color','k','linewidth',1)
            text(nn-7,-nn/7,['\theta = ' mint2str(theta(j)/pi*180) '°, \phi = ' mint2str(phi(j)/pi*180) '°'],'fontsize',3,'fontname','times');  
        end
        colormap(flipud(gray))
        print -dpng -r300 ModelPatterns
    end
    
    for j=1:size(mask,3) mask(:,:,j) = mask(:,:,j).*bck; mask(:,:,j) = mask(:,:,j)/sum(sum(mask(:,:,j))); end
    
    model = struct('rho',rho);
    model = setfield(model, 'theta', theta);
    model = setfield(model, 'phi', phi);
    model = setfield(model, 'mask', mask);
    model = setfield(model, 'fxx0', fxx0);
    model = setfield(model, 'fxx2', fxx2);
    model = setfield(model, 'fxz', fxz);
    model = setfield(model, 'byx0', byx0);
    model = setfield(model, 'byx2', byx2);
    model = setfield(model, 'byz', byz);
else
    nn = (size(model.mask,1)-1)/2;
    bck = Disk(nn); 
    rho = model.rho;
    theta = model.theta;
    phi = model.phi;
    mask = model.mask;
    fxx0 = model.fxx0;
    fxx2 = model.fxx2;
    fxz = model.fxz;
    byx0 = model.byx0;
    byx2 = model.byx2;
    byz = model.byz;
end

if nargin==0
    [file, pathname] = uigetfile('*.tif;*.tiff');
    filename{1} = [pathname file];
end

imnum = length(imfinfo(filename{1}));
tmp = mtiffread(filename{1});
for j=1:imnum
    im(:,:,j) = double(tmp(j).data);
end
clear tmp
if nargin<3 | isempty(roi)
    % localize region of interest
    clf; mim(sum(im,3));
    [b,a] = ginput(2);
    a=round(a); b=round(b);
    a(1) = max([a(1) 1]);
    a(end) = min([a(end) size(im,1)]);
    b(1) = max([b(1) 1]);
    b(end) = min([b(end) size(im,2)]);
    close
else
    a = roi(:,1);
    b = roi(:,2);
end

% find patterns
[x,y] = meshgrid(1:diff(b)+1,1:diff(a)+1);
[xx,yy] = meshgrid(-2:2,-2:2);
[xi,yi] = meshgrid(-2:0.05:2,-2:0.05:2);
[ind2,ind1] = meshgrid(1:size(xi,1),1:size(xi,2));
outfile = findstr(filename{1},'.');
outfile = filename{1}(1:outfile(end)-1);
tmp = dir([outfile '*.png']);
if length(tmp)>0
    outfile = [outfile mint2str(str2num(tmp(end).name(end-10:end-8))+1,3)];
else
    outfile = [outfile 'm001'];
end

clear tmp err sim xc yc z orient
close all

j = 0;
for k=1:length(filename)
    j0 = j;
    imnum = j0+length(imfinfo(filename{k}));
    for j=1:imnum
        [err, tmp, cim, sim] = FindPattern(im(a(1):a(2),b(1):b(2),j),mask,bck);
        tmp = cim./sqrt(err);
        tmp = tmp==max(tmp(:));
        xc = x(tmp); yc = y(tmp);
        xc = xc(end); yc = yc(end);
        if ~(b(1)-1+xc-nn-2>0 & b(1)-1+xc+nn+2<=size(im,2) & a(1)-1+yc-nn-2>0 & a(1)-1+yc+nn+2<=size(im,1))
            xc = []; yc =[];
            z(:,:,j0+j) = im(round(a(1)+diff(a)/2+(-nn:nn)),round(b(1)+diff(b)/2+(-nn:nn)),j);
            zz(:,:,j0+j) = 0*z(:,:,j0+j);
        end
        if ~isempty(xc)
            z(:,:,j0+j) = im(a(1)-1+(yc-nn:yc+nn),b(1)-1+(xc-nn:xc+nn),j);
            thetav = theta(sim(yc,xc))+(-2:2)*al_res/180*pi;
            phiv = phi(sim(yc,xc))+(-2:2)*be_res/180*pi;
            tmp = im(a(1)-1+(yc-nn-2:yc+nn+2),b(1)-1+(xc-nn-2:xc+nn+2),j);
            [tmporient, tmpzz, tmperr] = RefineOrientation(tmp,thetav,phiv,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
            zz(:,:,j0+j) = tmpzz(:,:,xx(tmperr==min(tmperr(:)))+3,yy(tmperr==min(tmperr(:)))+3);
            err = interp2(xx,yy,tmperr,xi,yi,'bicubic');
            i1 = ind1(err==min(err(:))); i1 = i1(1);
            i2 = ind2(err==min(err(:))); i2 = i2(1);
            coord(:,j0+j) = [b(1)-1+xc(end)+yi(i1,i2); a(1)-1+yc(end)+xi(i1,i2)];
            tmp = interp2(xx,yy,squeeze(tmporient(1,:,:)),xi,yi,'bicubic');
            orient(1,j0+j) = tmp(i1,i2);
            tmp = interp2(xx,yy,squeeze(tmporient(2,:,:)),xi,yi,'bicubic');
            orient(2,j0+j) = tmp(i1,i2);
        end
        j
    end
end

% format results
cnt = 1; bld = 1; res = [];
for j=1:size(z,3) 
    if sum(sum(zz(:,:,j)))>0
        res = [res; [j coord(:,j)' orient(1,j)/pi*180 orient(2,j)/pi*180]];
    end
    if pic
        subplot(4,5,cnt+5+5*floor((cnt-1)/5));
        mim(zz(:,:,j));
        if sum(sum(zz(:,:,j)))>0
            title({['\theta = ' mint2str(orient(1,j)/pi*180) '°']; ...
                ['\phi = ' mint2str(orient(2,j)/pi*180) '°']},'fontsize',12,'fontname','times');
        end
        subplot(4,5,cnt+5*floor((cnt-1)/5));
        mim(z(:,:,j));
        title(num2str(j),'fontsize',12);
        cnt = cnt+1;
        if cnt==11 | j==ceil(size(z,3))
            eval(['print -dpng ''' outfile '''-' mint2str(bld,3) '.png']);
            bld = bld+1;
            cnt = 1;
        end
    end
end

