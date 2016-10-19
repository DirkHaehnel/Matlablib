function [y, rad, yy, rr] = Circularize(tx, ty, z, nn)

% function Circularize calculates the distribution of a 2D field as a
% radial function from the position of its maximum

if isempty(tx) || isempty(ty)
    [tx,ty] = meshgrid(1:size(z,2),1:size(z,1));
end

if nargin<4 || isempty(nn)
    nn = 1e2;
end

ind = (z==max(z(:)));

if size(tx,1)==1 || size(tx,2)==1
    [tx,ty] = meshgrid(tx,ty);
end

rr = sqrt((tx-tx(ind)).^2 + (ty-ty(ind)).^2);
[rr ind] = sort(rr(:));
rad = (0:nn)/nn*max(rr(:));

y = mHist(rr(:),rad,z(:))./mHist(rr(:),rad);
yy = z(ind);

