function [xc, yc] = RefinePosition(errim,xc,yc,nn,precision)

if nargin<5 || isempty(precision)
    precision = 0.1;
end
v = -nn:nn;
vi = -nn:precision:nn;
[x,y] = meshgrid(v,v);
[xi,yi] = meshgrid(vi,vi);
        
for j=1:length(xc)
    erri = interp2(x,y,errim(xc(j)+v,yc(j)+v),xi,yi,'cubic');
    tmp = xi(erri==unique(min(erri(:))));
    xc(j) = xc(j) + tmp(1);
    tmp = yi(erri==unique(min(erri(:))));
    yc(j) = yc(j) + tmp(1);
end

