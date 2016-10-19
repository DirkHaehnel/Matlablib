function [res, z, a, b, mask] = Myosin(filename,mask,pic)

close all

NA = 1.45;
n0 = 1.52;
n = 1.333;
n1 = 1.333;
d0 = [];
d = 0; d1 = [];
lamex = 0.575;
mag = 160;
focus = 0.5;

be_res = 5; % minimum resolution of in-plane angle
al_res = 10; % minimum resolution of out-of-plane angle

if nargin<3 | isempty(pic)
    pic = true;
end

nn = 8;
n2 = 2*nn;
n3 = (2*nn+1)^2;
[x,y] = meshgrid(-nn:nn,-nn:nn);
p = angle(x+i*y);
r = sqrt(x.^2+y.^2);
bck = disk(nn); 

if nargin<2 | isempty(mask)
    [intx inty intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 3], 0, NA, n0, n, n1, d0, d, d1, lamex, mag, focus);
    rho = rho/16;
    cnt = 1;
    for k=0:al_res:90
        al = k/180*pi;
        if k==0
            jj = round(180/be_res);
            dbe = pi/jj;
        elseif k==90
            jj = 1;
            dbe = 0;
        else
            jj = round(cos(al)*360/be_res);
            dbe = 2*pi/jj;
        end
        for j=1:jj
            theta(cnt) = pi/2-al;
            be = dbe*(j-1);
            phi(cnt) = 2*pi-be;        
            int = real((cos(al)*(interp1(rho,fxx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,fxx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,fxz,r,'cubic')).*...
                conj(cos(al)*(interp1(rho,byx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,byx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,byz,r,'cubic')) + ...
                (cos(al)*sin(2*(p-be)).*interp1(rho,fxx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,fxz,r,'cubic')).*...
                conj(cos(al)*sin(2*(p-be)).*interp1(rho,byx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,byz,r,'cubic')));    
            mask(:,:,cnt) = int; cnt = cnt+1;
        end
        k
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
        print -dpng -r300 myosinmodel
    end
    
    for j=1:size(mask,3) mask(:,:,j) = mask(:,:,j).*bck; mask(:,:,j) = mask(:,:,j)/sum(sum(mask(:,:,j))); end
end

if nargin==0
    [filename, pathname] = uigetfile('*.tif;*.tiff');
    filename = [pathname filename]; 
end

% localize region of interest
imnum = length(imfinfo(filename));
tmp = mtiffread(filename);
for j=1:imnum
     im(:,:,j) = double(tmp(j).data); 
end
clear tmp
clf; mim(sum(im,3));
[b,a] = ginput(2);
a=round(a); b=round(b);
a(1) = max([a(1) 1]);
a(end) = min([a(end) size(im,1)]);
b(1) = max([b(1) 1]);
b(end) = min([b(end) size(im,2)]);
close

% find patterns
[xx,yy] = meshgrid(-2:2,-2:2);
[xi,yi] = meshgrid(-2:0.1:2,-2:0.1:2);
outfile = findstr(filename,'.');
outfile = filename(1:outfile(1)-1);
clear tmp err bim cim sim xc yc z orient intens
close all
for j=1:imnum
    [err, bim, cim, sim, xc, yc] = FindPattern(im(a(1):a(2),b(1):b(2),j),mask,bck,1,0.5,'mconv2(cim,disk(5))>tsh*sqrt(err)');
    if ~isempty(xc)
        tmp = interp2(xx,yy,err(yc-2:yc+2,xc-2:xc+2),xi,yi,'bicubic');
        coord(:,j) = [b(1)-1+xc(end)+xi(tmp==min(tmp(:))); a(1)-1+yc(end)+yi(tmp==min(tmp(:)))];
        z(:,:,j) = im(a(1)-1+(yc(end)-nn:yc(end)+nn),b(1)-1+(xc(end)-nn:xc(end)+nn),j);    
        [err, bla, c, s] = FindPattern(z(:,:,j),mask,bck);
        orient(j) = s((end+1)/2,(end+1)/2);
        intens(j) = c((end+1)/2,(end+1)/2);
    end
end

% format results
cnt = 1; bld = 1; res = [];
for j=1:size(z,3) 
    if orient(j)>0
        res = [res; [j coord(:,j)' theta(orient(j))/pi*180 phi(orient(j))/pi*180]];
        if pic
            subplot(4,5,cnt+5+5*floor((cnt-1)/5));                 
            mim(mask(:,:,orient(j))); 
            title({['\theta = ' mint2str(theta(orient(j))/pi*180) '°']; ...
                ['\phi = ' mint2str(phi(orient(j))/pi*180) '°']},'fontsize',12,'fontname','times');        
            subplot(4,5,cnt+5*floor((cnt-1)/5)); 
            mim(z(:,:,j)); 
            title(num2str(j),'fontsize',12); 
            cnt = cnt+1;
        end
    end;
    if pic
        if cnt==11 | j==ceil(size(z,3))
            eval(['print -djpeg ''' outfile '''-' mint2str(bld,2)]);
            bld = bld+1;
            cnt = 1;
        end
    end
end

