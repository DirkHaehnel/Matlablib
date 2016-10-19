function [err, z] = GaussRing(para,im)

[x,y] = meshgrid(1:size(im,2),1:size(im,1));

len = prod(size(im));

r = sqrt((x-para(1)).^2+(y-para(2)).^2);
%r = sqrt((x-para(1)).^2+para(6)*(y-para(2)).^2+para(7)*(x-para(1)).*(y-para(2)));

z = [ones(len,1) ...
        reshape(cumsum(ones(size(im))),len,1) ...
        reshape(cumsum(ones(size(im))')',len,1) ...
        reshape(exp(-r.^2/para(3)^2/2),len,1) ...
        reshape(exp(-(r-para(5)).^2/para(4)^2/2),len,1)];

z = reshape(z*lsqnonneg(z,reshape(im,prod(size(im)),1)),size(im,1),size(im,2));

mim(im,z); drawnow

err = sum(sum((im-z).^2));