% function [imm, xc, yc, ww] = Palm(im, ww)

close all 
fac = 30;
tsh = 1;

% read data
close all
[n,m,imlen] = size(data);

im = double(data);

% background reduction
[n,m,imlen] = size(data);
[x,y] = meshgrid(-(m-1)/2:(m-1)/2,-(n-1)/2:(n-1)/2);
bck = abs(ifft2(ifftshift(fftshift(fft2(mean(im,3))).*exp(-(x.^2+y.^2)/fac))));
im = im - repmat(bck,[1 1 imlen]);

% cut boundaries
boundary_width = 5;
a = round([boundary_width n-boundary_width]);
b = round([boundary_width m-boundary_width]);
im = im(a(1):a(2),b(1):b(2),:);
bck = bck(a(1):a(2),b(1):b(2));
[n,m,imlen] = size(im);

% determine size of single molecule image
% if nargin<2 || isempty(ww)
    mim(im(:,:,1))
    [b,a] = ginput(2);
    a = round(a);b=round(b);
    [mx,my,wx,wy] = GaussEllipse([],[],abs(im(a(1):a(2),b(1):b(2))-min(min(im(a(1):a(2),b(1):b(2))))));
    ww = sqrt(wx.*wy)/2;
% end

% find molecule positions
w3 = round(3*ww);
[xx,yy] = meshgrid(-w3:w3,-w3:w3);
mask = exp(-(xx.^2+yy.^2)/ww^2/2);
res = zeros(n,m);
for j=1:imlen
    [err, bim, cim, sim, xc, yc] = FindPattern(im(:,:,j)-min(min(im(:,:,j))), mask, [], [], 1, tsh);
    tmp = cim./sqrt(err).*(cim>tsh*sqrt(err));
    tmp(isnan(tmp)) = 0;
    res = res + tmp;
    mim(res)
end

