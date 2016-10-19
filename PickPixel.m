function [a,b] = PickPixel(im)

if nargin>0
    clf
    mim(im)
else
    figure(gcf)
end
[b,a] = ginput;
a = round(a);
b = round(b);

