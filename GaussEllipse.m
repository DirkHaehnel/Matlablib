function [mx, my, wx, wy, om, amp, zz] = GaussEllipse(tx, ty, z, pic)

if nargin==1
    z = tx;
end
if nargin==1 || isempty(tx)
    tx = (1:size(z,2))';
    ty = 1:size(z,1);
end
if (nargin>3 && ~isempty(pic)) || (nargin==2 && ~isempty(ty))
    pic = 1;
else
    pic = [];
end
tx = tx(:);
ty = ty(:)';
dtx = diff(tx(1:2));
dty = diff(ty(1:2));

if dtx<dty
    ti = min(ty):dtx:max(ty);
    z = interp1(ty,z,ti);
    ty = ti;
elseif dty<dtx
    ti = min(tx):dty:max(tx);
    z = interp1(tx,z',ti)';
    tx = ti(:);
end    

mx = sum(z.^2*tx)/sum(sum(z.^2));
my = sum(ty*z.^2)/sum(sum(z.^2));
wx = sqrt(sum(z*(tx-mx).^2)/sum(sum(z)))/2;
wy = sqrt(sum((ty-my).^2*z)/sum(sum(z)))/2;
if 1
    p = Simplex('Gauss',[mx wx],[-inf 0],[],[],[],tx,sum(z),1,[],pic);
    mx = p(1); wx = 2*p(2);
    p = Simplex('Gauss',[my wy],[-inf 0],[],[],[],ty',sum(z,2)',1,[],pic);
    my = p(1); wy = 2*p(2);
else
    opts = fitoptions('Method','Nonlinear','Normalize','On');
    g = fittype('c1+c2*exp(-2*(x-a)^2/b^2)','ind','x','prob',{'a','b'},'options',opts);
    set(opts,'Algorithm','Levenberg-Marquardt','Robust','on');
    tmp = fit(tx,sum(z)',g,'problem',{mx,wx});
    mx = tmp.a;
    wx = 4*tmp.b;
    tmp = fit(tx,sum(z)',g,'problem',{my,wy});
    my = tmp.a;
    wy = 4*tmp.b;
end

tmp = z;
[m,n] = size(tmp);
nn = floor(m*n/(m+1));
weight = [zeros(nn,n); ones(size(tmp))];
weight = reshape(weight,numel(weight),1);
weight = reshape(weight(1:(m+nn+1)*nn),m+nn+1,nn);
weight = sum(weight');
tmp = [zeros(nn,n); tmp];
tmp = reshape(tmp,numel(tmp),1);
tmp = reshape(tmp(1:(m+nn+1)*nn),m+nn+1,nn);
tmp = sum(tmp'); 
tmp(isnan(tmp)) = [];
md = sum((1:length(tmp)).*tmp)/sum(tmp);
p = sqrt(sum(((1:length(tmp))-md).^2.*tmp)/sum(tmp))/2;
if 1
    p = Simplex('Gauss',[md p],[-inf 0],[],[],[],[],tmp,1,weight,pic);
    wd = 2*min(dtx,dty)*p(2);
else
    % g = fittype({'c1+c2*exp(-2*(x-a)^2/b^2)+c3*weight'},'prob',{'a','b'},'options',opts);
    tmp = fit((1:length(tmp))',(tmp./(weight+(weight==0)))',g,'problem',{md,p});
    wd = 4*min(dtx,dty)*tmp.b;
end

om = -angle(1i*(wd.^2-wx.^2-wy.^2)+(wx.^2-wy.^2))/2;
v = sqrt([cos(om)^2 sin(om)^2; sin(om)^2 cos(om)^2; 1-sin(2*om) 1+sin(2*om)]\[wx^2;wy^2;wd^2]);

wx = v(1);
wy = v(2);
if om>pi
    om=om-pi;
end

[xx,yy] = meshgrid(tx,ty);
zz = exp(-2*((xx-mx)*cos(om)+(yy-my)*sin(om)).^2/wx^2-2*(-(xx-mx)*sin(om)+(yy-my)*cos(om)).^2/wy^2);
amp = [ones(numel(xx),1) zz(:)]\z(:);
zz = reshape([ones(numel(xx),1) zz(:)]*amp,size(z,1),size(z,2));
amp = amp(2)*wx*wy/2;

phi = 0:pi/100:2*pi;
pcolor(tx,ty,z); axis image; shading flat
hold on
plot(mx+diff(tx(1:2))/2+wx*cos(phi)*cos(om)-wy*sin(phi)*sin(om),my+diff(ty(1:2))/2+wx*cos(phi)*sin(om)+wy*sin(phi)*cos(om),'y','linewidth',1);
hold off

