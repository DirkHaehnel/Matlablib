function [err, z, zz] = Gauss2DFun(p, tx, ty, z, pic)

if nargin==2
    z = tx;
end
if nargin==2 || isempty(tx)
    [tx,ty] = meshgrid(1:size(z,2),1:size(z,1));
end
if (nargin>4 && ~isempty(pic)) 
    pic = 1;
else
    pic = [];
end
if length(p)==3
    p(4) = p(3);
end
z(~isfinite(z)) = 0;

zz = exp(-(tx-p(1)).^2/2/p(3)^2-(ty-p(2)).^2/2/p(4)^2);

c = [ones(numel(zz),1) zz(:)]\z(:);
zz = c(1)*ones(size(z)) + c(2)*zz;

CombineImages(cat(3,z,zz),1,2);

err = sum((z(:)-zz(:)).^2./abs(zz(:)))

