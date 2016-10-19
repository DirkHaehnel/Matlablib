function [im, xc, yc, amp, im0, bck] = FindMolecules(name)

if nargin==0
    [filename, pathname] = uigetfile('*.tif', 'Interactive mode data:', 0, 0);
    name = [pathname filename];
end

fac = 30;
mol_width = 2;
tsh = 0.3;

close all
imlen = length(imfinfo(name));
for j=1:imlen
    tmp = imread(name,j);
    im(:,:,j) = double(tmp);
end
if strcmp(name(end-3:end),'.tif')
    name = name(1:end-4);
end
    
clear tmp
[n,m,imlen] = size(im);
a = round([2 n-2]);
b = round([2 m-2]);
im = im(a(1):a(2),b(1):b(2),:);
im0 = im;

[x,y] = meshgrid(-size(im,2)/2+0.5:size(im,2)/2,-size(im,1)/2+0.5:size(im,1)/2);
bck = abs(ifft2(ifftshift(fftshift(fft2(sum(im,3))).*exp(-(x.^2+y.^2)/fac))))/imlen;
for j=1:imlen
    im(:,:,j) = im(:,:,j)-bck;
end

ww = 3*ceil(mol_width);
[x,y] = meshgrid(-ww:ww,-ww:ww);
mask = exp(-(x.^2+y.^2)/mol_width^2);
[tmp, tmp, tmp, tmp, xc, yc] = FindPattern(sum(im,3)-min(min(sum(im,3))), mask, ones(size(x)), 1, tsh); 
for j=1:imlen
    [tmp, tmp, cim(:,:,j)] = FindPattern(im(:,:,j)-min(im(:)), mask, ones(size(x)), 1, tsh); 
    for k=1:length(xc)
        amp(j,k) = cim(yc(k),xc(k),j);
    end
end

