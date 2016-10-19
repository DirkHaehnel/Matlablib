function [err, c, z] = Gauss2fFIMDAFull(para,velo,w0,a0,av,lam,delta,md,tt,yd)

% para(1:end/2) - concentration
% para(end/2+1:end) - brigthness

NA = 6.0221e+023;

if mod(length(para),2)==0
    fac = 1;
else
    fac = para(end);
    para(end) = [];
end
bg = para(end-1:end);
para = para(1:end-2);
cc = para(1:end/2)*NA/1e24;
bb = para(end/2+1:end);
len = length(para)/2;

xi = exp(i*(0:max(md))/(max(md)+1)*2*pi)-1;

Nz = 2; %1e1;
zmax = 1e1*a0; %evaluate function until amplitude function has fallen off to 1%
z = (0.5:Nz)/Nz*zmax;
Nx = 1e2;
rho = (0.5:Nx)/Nx*4;
[x,y] = meshgrid(rho,rho);
ww = LaserBeamFit([w0 0],lam(1),z); 
weight = ww.^2/2*4*zmax/Nx/Nz;
s2 = sqrt(2);
kappa = sqrt(pi/2)/2*w0^2*D1CEF([a0 0],av,lam(2),z)./ww/velo;

for j1=1:length(xi)
    for j2=1:length(xi)
        tmp = bg(1)*xi(j1) + bg(2)*xi(j2);
        for jz=1:length(z)
            uu = kappa(jz)*(erf(s2*(velo*tt-rho)/ww(jz))-erf(s2*rho/ww(jz)));
            u1 = exp(-(rho-delta/2)'.^2/ww(jz)^2);
            u2 = exp(-(rho+delta/2)'.^2/ww(jz)^2);
            for k=1:len
                tmp = tmp - cc(k)*sum(sum(1 - exp(bb(k)*(xi(j1)*u1*uu+fac*xi(j2)*u2*uu))))*weight(jz);
            end
        end
        g(j1,j2) = tmp;
    end
    j1
end
g = exp(g);
g = abs(real(fft2(g)));

if nargin>9
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

% [g, md] = Gauss2fFIMDAFull([1e-9 300 1 1],1e3,350,150,100e3,[640 670]/1.33,400,0:100,3e-3); mim(g)