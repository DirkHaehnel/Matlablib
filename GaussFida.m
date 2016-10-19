function [err, c, z] = GaussFida(para,w0,a0,av,lam,md,yd)

% p(1:end/2) - concentration
% p(end/2+1:end) - brigthness

NA = 6.0221e+023;

if mod(length(para),2)==1
    bg = para(end);
    para = para(1:end-1);
else
    bg = 0;
end
para(1:end/2) = para(1:end/2)*NA/1e24;

xi = exp(i*(0:max(md))/(max(md)+1)*2*pi)-1;
len = length(para);

Nz = 1e3;
zmax = 1e2*a0; %evaluate function until amplitude function has fallen off to 1%
z = (0.5:Nz)/Nz*zmax;
Nx = 1e2;
rho = (0.5:Nx)/Nx*4;
w2 = LaserBeamFit([w0 0],lam(1),z).^2;
uu = w0^2*exp(-rho'.^2)*(D1CEF([a0 0],av,lam(2),z)./w2);
weight = 4*pi*rho'*w2/2*4*zmax/Nx/Nz;

for j=1:length(xi)
    tmp = bg*xi(j);
    for k=1:len/2
       tmp = tmp - para(k)*sum(sum((1 - exp(xi(j)*para(end/2+k)*uu)).*weight));
    end
    g(j) = exp(tmp);
end
g = abs(real(fft(g)));

if nargin>6
    md = md(:); yd = yd(:); g = g(:);
    c = lsqnonneg(g,yd);
    z = c*g;
    loglog(md,yd,'o',md,z);
    drawnow
    err = sum((yd(yd>0)-z(yd>0)).^2./abs(yd(yd>0)));
else
    err = g;
    c = 0:md;
end

