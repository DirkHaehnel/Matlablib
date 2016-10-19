function [im, cim, resx, resy, im0, bck, err] = Iris(name, mol_width);

% converting image stacks into movies

if nargin==0
    [filename, pathname] = uigetfile('*.tif', 'Interactive mode data:', 0, 0)
    name = [pathname filename];
end

fac = 30;
if nargin<2 || isempty(mol_width)
    mol_width = 2;
end
tsh = 0.5;

close all
tmp = tiffread(name);
for j=1:length(tmp)
    im(:,:,j) = double(tmp(j).data);
end
if strcmp(name(end-3:end),'.tif')
    name = name(1:end-4);
end
    
clear tmp
% figure
% set(gcf,'units','normalized','position',[0 0 0.9 0.9])
% mim(log(tmp));
% [b,a] = ginput(2);
% a = round(a); b = round(b);
[n,m,imlen] = size(im);
% a = round([0.1 0.9]*n);
% b = round([0.1 0.9]*m);
% im = im(a(1):a(2),b(1):b(2),:);
im0 = im;

[x,y] = meshgrid(-size(im,2)/2+0.5:size(im,2)/2,-size(im,1)/2+0.5:size(im,1)/2);
bck = abs(ifft2(ifftshift(fftshift(fft2(sum(im,3))).*exp(-(x.^2+y.^2)/fac))))/imlen;
for j=1:imlen
    im(:,:,j) = im(:,:,j)-bck;
end

[x,y] = meshgrid(-3*ceil(mol_width):3*ceil(mol_width),-3*ceil(mol_width):3*ceil(mol_width));
mask = exp(-(x.^2+y.^2)/mol_width^2);
for j=1:imlen
    [err(:,:,j), tmp, cim(:,:,j), tmp, xc, yc] = FindPattern(im(:,:,j)-min(im(:)), mask, ones(size(x)), 1, tsh); 
    resx{j}=xc; 
    resy{j}=yc; 
end

return

imlen = size(im,3);

cx = [min(im(:)) max(im(:))];
for j=1:imlen
    mim(im(:,:,j));
    caxis(cx);
    mov(j) = getframe;
end
movie2avi(mov,[name '_raw.avi'],'compression','none');

cx = [min(im(:)) max(im(:))];
for j=1:imlen
    mim(cim(:,:,j));
    caxis(cx);
    hold on;
    plot(resx{j},resy{j},'ob');
    hold off
    mov(j) = getframe;
end
movie2avi(mov,[name '.avi'],'compression','none');

bin = 0:ceil(size(im,1)^2+size(im,2).^2);
dist = zeros(length(bin),imlen-1);
for k=1:imlen-1
    cnt = 0;
    for j=1:imlen-k
        tmp = (resx{j}'*ones(1,length(resx{j+k}))-ones(length(resx{j}),1)*resx{j+k}).^2 + ...
            (resy{j}'*ones(1,length(resy{j+k}))-ones(length(resy{j}),1)*resy{j+k}).^2;
        dist(:,k) = dist(:,k) + mHist(tmp(:),bin);
        cnt = cnt + prod(size(tmp));
    end
    dist(:,k) = dist(:,k)/cnt;
end
c = polyfit(1:imlen-1,bin*dist./sum(dist),1);
plot(1:imlen-1,bin*dist./sum(dist),1:imlen-1,polyval(c,1:imlen-1));
xlabel('frame #');
ylabel('\langle\itr\rm^2\rangle [pix^2]')
title(['slope = ' mnum2str(c(1),2,1) ' pix^2/T_{frame}'])
eval(['print -dpng ''' name '''.png' ])

