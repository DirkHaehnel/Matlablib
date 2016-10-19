function err = EllipseFit(para,im)

[m,n] = size(im);
[x,y] = meshgrid(1:n,1:m);

z = sqrt((cos(para(3))*(x-para(1))-sin(para(3))*(y-para(2))).^2+(sin(para(3))*(x-para(1))+cos(para(3))*(y-para(2))).^2*para(4));
z = z>=para(5) & z<=para(6);
col = ones(prod(size(im)),1);
c = [col z(:)]\im(:);

mim(z.*im); drawnow;

err = sum((im(:) - [col z(:)]*c).^2);
