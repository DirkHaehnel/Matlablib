function [res, z, zz, a, b, model] = PatternMatching(file,pic,roi,model)

close all

NA = 1.49; 		% Numerical aperture of objective 
n0 = 1.52;		% refractive index of immersion medium
n = 1.0;		% refractive index of sample
n1 = 1.0;		% refractive index of topping medium
d0 = [];				
d = 0.01; 
d1 = [];
lamex = 0.57;   % excitation wavelength in mum
mag = 100;		% magnification
focus = 0.7;		% defcousing in mum
pixel = 16;		% CCD pixel size in mum

be_res = 10; % minimum resolution of in-plane angle
al_res = 10; % minimum resolution of out-of-plane angle

if nargin==0
    [file, pathname] = uigetfile('*.tif;*.tiff');
    filename{1} = [pathname file];
elseif ~iscell(file)
    filename{1} = file;
else
    filename = file;
end

if nargin<2 || isempty(pic)
    pic = true;
end

if nargin<4 || isempty(model)
    nn = 15;
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
    
    [intx inty intz, rho, tmp, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 3], 0, NA, n0, n, n1, d0, d, d1, lamex, mag, focus);
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

if ~file==0
    if ~isnumeric(filename)
        imnum = length(imfinfo(filename{1}));
        if nargin<3 || isempty(roi)
            % localize region of interest
            if imnum==1
                im = double(imread(filename{1}));
            else
                im = double(imread(filename{1},1));
                for j=2:imnum
                    im = im + double(imread(filename{1},j));
                end
            end
            clear tmp
            clf; mim(im);
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
        [x,y] = meshgrid(1:diff(a)+1,1:diff(b)+1);
        [xx,yy] = meshgrid(-2:2,-2:2);
        [xi,yi] = meshgrid(-2:0.05:2,-2:0.05:2);
        [ind2,ind1] = meshgrid(1:size(xi,1),1:size(xi,2));
        outfile = strfind(filename{1},'.');
        outfile = filename{1}(1:outfile(end)-1);
        tmp = dir([outfile '*.png']);
        if ~isempty(tmp)
            outfile = [outfile mint2str(str2num(tmp(end).name(end-10:end-8))+1,3)];
        else
            outfile = [outfile 'm001'];
        end
        clear tmp
    else
        pic = false;
    end
    
    close all
    
    j = 0;
    if isnumeric(filename)
        kmax = 1;
    else
        kmax = ength(filename);
    end
    for k=1:kmax
        j0 = j;
        if isnumeric(filename)
            imnum = size(filename,3);
        else
            imnum = length(imfinfo(filename{k}));
        end
        for j=1:imnum
            if isnumeric(filename)
                im = filename(:,:,k);
            else
                if imnum==1
                    im = double(imread(filename{k}));
                else
                    im = double(imread(filename{k},j));
                end
            end
            % [tmp, tmp, tmp, sim, xc, yc] = FindPattern(im(a(1):a(2),b(1):b(2)),mask,bck,bck,1,1,'mconv2(cim,Disk(5))>tsh*sqrt(err)');
            
            a(1) = max([a(1) 3]);
            a(2) = min([a(2) size(im,1)-2]);
            b(1) = max([b(1) 3]);
            b(2) = min([b(2) size(im,2)-2]);
            nn   = min([floor(diff(a)/2) floor(diff(b)/2)]);
            
            [err, tmp, cim, sim, xc, yc] = FindPattern(im(a(1):a(2),b(1):b(2)),mask,bck,bck);
            if ~isempty(xc)
                if abs(xc-diff(b)/2-1)>4 || abs(yc-diff(a)/2-1)>4
                    xc = []; yc  = [];
                        z(:,:,j0+j) = im(a(1)+round(diff(a)/2)+(-nn:nn),b(1)+round(diff(b)/2)+(-nn:nn));
                    zz(:,:,j0+j) = 0*z(:,:,j0+j);
                end
            end
            if ~isempty(xc)
                z(:,:,j0+j) = im(a(1)-1+(yc(end)-nn:yc(end)+nn),b(1)-1+(xc(end)-nn:xc(end)+nn));
                thetav = theta(sim(yc(end),xc(end)))+(-2:2)*al_res/180*pi;
                phiv = phi(sim(yc(end),xc(end)))+(-2:2)*be_res/180*pi;
                tmp = im(a(1)+(yc(end)-nn-2:yc(end)+nn+2),b(1)+(xc(end)-nn-2:xc(end)+nn+2));
                [tmporient, tmpzz, tmperr] = RefineOrientation(tmp,thetav,phiv,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
                zz(:,:,j0+j) = tmpzz(:,:,yy(tmperr==min(tmperr(:)))+3,xx(tmperr==min(tmperr(:)))+3);
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
            if (cnt==11 || j==ceil(size(z,3))) 
                eval(['print -dpng ''' outfile '''-' mint2str(bld,3) '.png']);
                bld = bld+1;
                cnt = 1;
            end
        end
    end
else
   res = []; 
   z = [];
   zz = [];
   a = [];
   b = []; 
end

return

% anything below this line is not executed and serves only as repository of ideas
plot(res(:,1),res(:,5),res(:,1),unwrap(res(:,5)/180*2*pi)/2/pi*180)
plot(res(:,1),90-abs(90-res(:,4)),res(:,1),unwrap(res(:,5)/180*2*pi)/2/pi*180)

theta=(90-abs(90-res(:,4)))/180*pi; phi=unwrap(res(:,5)/180*2*pi)/2;
clf
sphere(100); shading interp; 
h=line([0 0],[0 0],[0 1.]); set(h,'color','c','linewidth',1); 
h=line([0 1.],[0 0],[0 0]); set(h,'color','c','linewidth',1); 
h=line([0 0],[0 0],[1.0 1.2]); set(h,'color','b','linewidth',1); 
h=line([1.0 1.2],[0 0],[0 0]); set(h,'color','b','linewidth',1); 
alpha(0.2); camlight(120, 60); axis image; axis off;
hold on
plot3(sin(theta).*cos(phi),sin(theta).*sin(phi),cos(theta),'linewidth',1);
hold off
view([120 30])
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2])
axis vis3d
cameratoolbar


theta=(90-abs(90-res(:,4)))/180*pi; phi=unwrap(res(:,5)/180*2*pi)/2;
t = 1:size(res,1);
t = t(~(res([1:end-1],1)+1==res([2:end],1)));
t = [0 t size(res,1)];
for j=1:length(t)-1
    pp(j).theta = theta(t(j)+1:t(j+1));
    pp(j).phi = phi(t(j)+1:t(j+1));
end
clf
sphere(100); shading interp; 
h=line([0 0],[0 0],[0 1.]); set(h,'color','c','linewidth',1); 
h=line([0 1.],[0 0],[0 0]); set(h,'color','c','linewidth',1); 
h=line([0 0],[0 0],[1.0 1.2]); set(h,'color','b','linewidth',1); 
h=line([1.0 1.2],[0 0],[0 0]); set(h,'color','b','linewidth',1); 
alpha(0.2); camlight(120, 60); axis image; axis off;
hold on
for j=1:length(pp)
    plot3(sin(pp(j).theta).*cos(pp(j).phi),sin(pp(j).theta).*sin(pp(j).phi),cos(pp(j).theta),'linewidth',1,'color',[j/length(pp) 0 1-j/length(pp)]);
end
hold off
view([120 30])
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2])
axis vis3d
cameratoolbar


theta=(90-abs(90-res(:,4)))/180*pi; phi=unwrap(res(:,5)/180*2*pi)/2;
tmp = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
bin = min(sum(tmp(1:end-1,:).*tmp(2:end,:),2)):0.001:1;
j = 1;
tmp = [sin(pp(j).theta).*cos(pp(j).phi), sin(pp(j).theta).*sin(pp(j).phi), cos(pp(j).theta)];
hh = mhist(sum(tmp(1:end-1,:).*tmp(2:end,:),2),bin);
for j=2:length(pp)
    tmp = [sin(pp(j).theta).*cos(pp(j).phi), sin(pp(j).theta).*sin(pp(j).phi), cos(pp(j).theta)];
    hh = hh + mhist(sum(tmp(1:end-1,:).*tmp(2:end,:),2),bin);
end
plot(acos(bin)/pi*180,hh)
