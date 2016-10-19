function [err, c, z] = GaussFida(para,w0,a0,av,lam,md,yd)

% p(1:end/2) - concentration
% p(end/2+1:end) - brigthness

if mod(length(para),2)==1
    bg = para(end);
    para = para(1:end-1);
else
    bg = 0;
end

xi = exp(i*(0:max(md)-1)/max(md)*2*pi)-1;
len = length(para);

Nz = 1e3;
zmax = 1e2*a0; %evaluate function until amplitude function has fallen off to 1%
z = (0.5:Nz)/Nz*zmax;
Nx = 1e2;
rho = (0.5:Nx)/Nx*4;
w2 = LaserBeamFit([w0 0],lam(1),z).^2;
uu = w0^2*exp(-rho'.^2)*(D1CEF([a0 0],av,lam(2),z)./w2);
weight = 4*pi*rho'*w2/2*4*zmax/Nx/Nz;

md = md(:)';
mi = 0.1:0.1:max(md); len = max(md)+1;
mm = 0:max(md); m1 = mm; m1(1)=1;
g = 0;
for k=1:length(para)/2
    g0 = mhist(para(end/2+k)*uu(:),mi,weight(:));
    mm = 0:max(md); m1 = mm; m1(1)=1;
	g0 = g0'*(exp(log(mi')*mm-mi'*ones(1,length(mm))-ones(length(mi),1)*cumsum(log(m1))));
    fac = para(k);
    g1 = fac*g0;
    g0 = flipud(g0);
    for j=2:10
        tmp = conv(g0,g1);
        fac = para(k)*fac/j;
        g1 = g1 + fac*tmp(1:length(g0));
    end
    g = g + g1(md+1);
end
   
if nargin>6
    md = md(:); yd = yd(:); g = g(:);
    c = lsqnonneg(g,yd);
    z = c*g;
    loglog(md,yd,'o',md,z);
    drawnow
    err = sum((yd(yd>0)-z(yd>0)).^2./abs(yd(yd>0)));
else
    err = g;
    c = md;
end

