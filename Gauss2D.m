function [mx, my, wx, wy, amp, xx, yy, zz] = Gauss2D(tx, ty, z, om, pic)

if nargin==1
    z = tx;
end
if nargin==1 || isempty(tx)
    tx = (1:size(z,2))';
    ty = 1:size(z,1);
end
if (nargin>4 && ~isempty(pic)) 
    pic = 1;
else
    pic = [];
end
tx = tx(:);
ty = ty(:)';
dtx = diff(tx(1:2));
dty = diff(ty(1:2));

zz = z; zz(~isfinite(z)) = 0;
mx = sum(zz.^2*tx)/sum(sum(zz.^2));
my = sum(ty*zz.^2)/sum(sum(zz.^2));
wx = sqrt(sum(zz*(tx-mx).^2)/sum(sum(zz)))/2;
wy = sqrt(sum((ty-my).^2*zz)/sum(sum(zz)))/2;
p = Simplex('Gauss',[mx wx],[-inf 0],[],[],[],tx,sum(z),1,[],pic);
mx = p(1); wx = 2*p(2);
p = Simplex('Gauss',[my wy],[-inf 0],[],[],[],ty',sum(z,2),1,[],pic);
my = p(1); wy = 2*p(2);

if nargin>3 && ~(isempty(om) || om==0)
    tmp = z;
    [m,n] = size(tmp);
    nn = floor(m*n/(m+1));
    weight = [zeros(nn,n); ones(size(tmp))];
    weight = reshape(weight,numel(weight),1);
    weight = reshape(weight(1:(m+nn+1)*nn),m+nn+1,nn);
    weight = sum(weight,2);
    tmp = [zeros(nn,n); tmp];
    tmp = reshape(tmp,numel(tmp),1);
    tmp = reshape(tmp(1:(m+nn+1)*nn),m+nn+1,nn);
    tmp = sum(tmp,2);
    tmp(isnan(tmp)) = [];
    md = sum((1:length(tmp))'.*tmp)/sum(tmp);
    p = sqrt(sum(((1:length(tmp))'-md).^2.*tmp)/sum(tmp))/2;
    p = Simplex('Gauss',[md p],[-inf 0],[],[],[],[],tmp,1,weight,pic);
    wd = 2*min(dtx,dty)*p(2);

    v = sqrt([cos(om)^2 sin(om)^2; sin(om)^2 cos(om)^2; 1-sin(2*om) 1+sin(2*om)]\[wx^2;wy^2;wd^2]);
    wx = v(1);
    wy = v(2);
else
    om = 0;
end

[xx,yy] = meshgrid(tx,ty);
zz = exp(-2*(xx-mx).^2/wx^2-2*(yy-my).^2/wy^2);
col = ones(numel(xx),1);
ind = isfinite(z(:));
amp = [col(ind) zz(ind)]\z(ind);
amp(2) = amp(2)*wx*wy/2;

phi = 0:pi/100:2*pi;
pcolor(tx,ty,z); axis image; shading flat
hold on
plot(mx+diff(tx(1:2))/2+wx*cos(phi)*cos(om)-wy*sin(phi)*sin(om),my+diff(ty(1:2))/2+wx*cos(phi)*sin(om)+wy*sin(phi)*cos(om),'y','linewidth',1);
hold off

