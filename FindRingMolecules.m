function [im, xc, yc, sc, mask, err, cim, sim] = FindRingMolecules(name,flag);

if nargin==0
    [name, path] = uigetfile('*.tif');
    name = [path name];
end

if isstr(name)
    imlen = length(imfinfo(name));
    im = double(imread(name,1)); 
    mm = [min(im(:)) max(im(:))];
    for j=2:imlen 
        tmp = double(imread(name,j));
        mm(1) = min([mm(1) min(tmp(:))]);
        mm(2) = max([mm(2) max(tmp(:))]); 
        im = im + tmp;
    end
else
    im = name;
end

[int, int, int, rho, int, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 3], 0, 1.2, 1.51, 1, 1, [], 0, [], 0.65, 60, 1.0);
rho = rho/6.45;
nn = 10;
n2 = 2*nn;
n3 = (2*nn+1)^2;
[x,y] = meshgrid(-nn:nn,-nn:nn);
p = angle(x+i*y);
r = sqrt(x.^2+y.^2);

cnt = 1;
for j=0:3
    al = j/3*pi/2;
    int = real((cos(al)*(interp1(rho,fxx0,r,'cubic')+cos(2*p).*interp1(rho,fxx2,r,'cubic'))+sin(al)*cos(p).*interp1(rho,fxz,r,'cubic')).*...
        conj(cos(al)*(interp1(rho,byx0,r,'cubic')+cos(2*p).*interp1(rho,byx2,r,'cubic'))+sin(al)*cos(p).*interp1(rho,byz,r,'cubic')) + ...
        (cos(al)*sin(2*p).*interp1(rho,fxx2,r,'cubic')+sin(al)*sin(p).*interp1(rho,fxz,r,'cubic')).*...
        conj(cos(al)*sin(2*p).*interp1(rho,byx2,r,'cubic')+sin(al)*sin(p).*interp1(rho,byz,r,'cubic')));    
    switch j
        case 0,
            al = 0;
            kmax = 15;
            dbe = pi/15;
        case 1,
            al = 30/180*pi;
            kmax = 12;
            dbe = pi/6;
        case 2,
            al = 60/180*pi;
            kmax = 8;
            dbe = pi/4;
        case 3
            al = pi/2;
            kmax = 1;
            dbe = 0;
    end
    for k = 1:kmax 
        theta(cnt) = pi/2-al;
        be = dbe*(k-1)/pi*180;
        phi(cnt) = be;        
        mask(:,:,cnt) = imrotate(int,be,'bicubic','crop');
        cnt = cnt+1        
    end
end
for j=1:size(mask,3) subplot(5,9,9+j); mim(mask(:,:,j)); title(num2str(j)); end

bck = disk(nn); 
for j=1:size(mask,3) mask(:,:,j) = mask(:,:,j).*bck; end
[err, cim, sim, xc, yc, sc] = FindPattern(im,mask,bck,11);

if nargin==2 & flag>0
    save(name(1:end-4),'im','err','cim','sim','xc','yc','sc','mask');
end

break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf; for j=1:size(mask,3) subplot(5,9,9+j); mim(mask(:,:,j)); title(num2str(j)); end
for j=1:length(xc)
    subplot(5,9,1); 
    tmp=im(yc(j)+(-nn:nn),xc(j)+(-nn:nn)); 
    mim(tmp); 
    title([num2str(j) '  ' num2str(sqrt(err(yc(j),xc(j)))/cim(yc(j),xc(j)))]);
    subplot(5,9,9+sc(j)); 
    hold 
    plot([1 2*nn+1.5 2*nn+1.5 1 1],[1 1 2*nn+1.5 2*nn+1.5 1],'y','linewidth',2);
    hold
    pause; 
    mim(mask(:,:,sc(j))); title(num2str(sc(j)));
end

