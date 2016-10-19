function [a,b]=fftdist(im)

im = ifftshift(abs(ifft2(im)));
[x,y] = meshgrid(1:size(im,2),1:size(im,1));
x0 = x(im==max(max(im)));
y0 = y(im==max(max(im)));
r = sqrt((x-x0).^2+(y-y0).^2);
[a,b] = mHist(r(:),0:max(r(:)),im(:));
b(1) = [];
a = a(2:end)./b';

