function [err, c, z] = Gauss2fFIMDAMCS(para,velo,w0,a0,av,lam,delta)

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
[xi,eta] = meshgrid(xi,xi);

F = @(x) D1CEF(a0,av,lam(2),x);
v0 = pi/2*w0^2*quadl(F,0,1e10);

F = @(x) D1CEF(a0,av,lam(2),x).^2.*exp(-delta^2./LaserBeamFit(w0,lam(1),x).^2).*...
    (sqrt(pi)*velo*tt./LaserBeamFit(w0,lam(1),x).*erf(velo*tt./LaserBeamFit(w0,lam(1),x))+...
    exp(-(velo*tt./LaserBeamFit(w0,lam(1),x)).^2)-1);
g12 = pi/8/velo^2*w0^4*quadl(F,0,1e10);

F = @(x) D1CEF(a0,av,lam(2),x).^2.*...
    (sqrt(pi)*velo*tt./LaserBeamFit(w0,lam(1),x).*erf(velo*tt./LaserBeamFit(w0,lam(1),x))+...
    exp(-(velo*tt./LaserBeamFit(w0,lam(1),x)).^2)-1);
g11 = pi/8/velo^2*w0^4*quadl(F,0,1e10);

g = tt*(bg(1)*xi + bg(2)*eta);
for j=1:len
    g = g + cc(j)*tt*v0*bb(j)*(xi + fac*eta);
    g = g + cc(j)*bb(j)^2*(g11*(xi.^2 + fac^2*eta.^2) + 2*g12*fac*xi.*eta);
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

