function err = RingFit(para,im)

[m,n] = size(im);
[x,y] = meshgrid(1:n,1:m);

z = sqrt((x-para(1)).^2+(y-para(2)).^2);
z = z>=para(3) & z<=para(3)+para(4);
col = ones(numel(im),1);
c = [col z(:)]\im(:);

mim(z.*im); drawnow;

err = sum((im(:) - [col z(:)]*c).^2);
