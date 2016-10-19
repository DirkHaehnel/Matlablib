function [a,b] = PickVec(im)

if nargin>0
    clf
    mim(im)
else
    figure(gcf)
end
[b,a] = ginput(2);
a = round(a);
b = round(b);
a = a(1):a(2);
b = b(1):b(2);

