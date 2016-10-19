% calculate the FRET Imaging of single molecules with NV centers

lambda = 750;
zd = 0;
za = 2.5;
rhoa = (0.1:0.2:20);
n0 = 2.42;
nd = n0;
ni = 1;
na = 1.46;
n1 = na;
d0 = [];
dd = 5;
di = 2;
da = 5;
d1 = [];

[ex,ez,ex0,ez0] = DipoleField(zd/lambda*2*pi,rhoa/lambda*2*pi,za/lambda*2*pi,n0,nd,ni,na,n1,d0,dd/lambda*2*pi,di/lambda*2*pi,da/lambda*2*pi,d1);

phi = (0:360)/180*pi;
one = ones(size(phi));

% fields from x-dipole
xx = ex(1,:).'*one + ex(2,:).'*cos(2*phi);
xy = ex(2,:).'*sin(2*phi);
xz = ex(3,:).'*cos(phi);

% fields from y-dipole
yy = ex(1,:).'*one - ex(2,:).'*cos(2*phi);
yx = -ex(2,:).'*sin(2*phi);
yz = ex(3,:).'*sin(phi);

% fields from z-dipole
zx = ez(2,:).'*cos(phi);
zy = ez(2,:).'*sin(phi);
zz = ez(1,:).'*one;

mm = max(max(abs(yx).^2+abs(zx).^2));
mm = max([mm max(max(abs(yy).^2+abs(zy).^2))]);
mm = max([mm max(max(abs(yz).^2+abs(zz).^2))]);
mm = max([mm max(max(abs(xx).^2+abs(yx).^2))]);
mm = max([mm max(max(abs(xy).^2+abs(yy).^2))]);
mm = max([mm max(max(abs(xz).^2+abs(yz).^2))]);

% x-acceptor & x-dark axis donor
subplot(2,3,1)
pcolor(rhoa'*sin(phi),rhoa'*cos(phi),abs(yx).^2+abs(zx).^2); axis image; shading interp; set(gca,'fontsize',10)
xlabel('\itx\rm (nm)','fontsize',10)
ylabel('\ity\rm (nm)','fontsize',10)
caxis([0 mm]);
title('x-acceptor & x-dark axis donor','fontsize',10)
% y-acceptor & x-dark axis donor
subplot(2,3,2)
pcolor(rhoa'*sin(phi),rhoa'*cos(phi),abs(yy).^2+abs(zy).^2); axis image; shading interp; set(gca,'fontsize',10)
xlabel('\itx\rm (nm)','fontsize',10)
ylabel('\ity\rm (nm)','fontsize',10)
caxis([0 mm]);
title('y-acceptor & x-dark axis donor','fontsize',10);
% z-acceptor & x-dark axis donor
subplot(2,3,3)
pcolor(rhoa'*sin(phi),rhoa'*cos(phi),abs(yz).^2+abs(zz).^2); axis image; shading interp; set(gca,'fontsize',10)
xlabel('\itx\rm (nm)','fontsize',10)
ylabel('\ity\rm (nm)','fontsize',10)
caxis([0 mm]);
title('z-acceptor & x-dark axis donor','fontsize',10);
% x-acceptor & z-dark axis donor
subplot(2,3,4)
pcolor(rhoa'*sin(phi),rhoa'*cos(phi),abs(xx).^2+abs(yx).^2); axis image; shading interp; set(gca,'fontsize',10)
xlabel('\itx\rm (nm)','fontsize',10)
ylabel('\ity\rm (nm)','fontsize',10)
caxis([0 mm]);
title('x-acceptor & z-dark axis donor','fontsize',10);
% y-acceptor & z-dark axis donor
subplot(2,3,5)
pcolor(rhoa'*sin(phi),rhoa'*cos(phi),abs(xy).^2+abs(yy).^2); axis image; shading interp; set(gca,'fontsize',10)
xlabel('\itx\rm (nm)','fontsize',10)
ylabel('\ity\rm (nm)','fontsize',10)
caxis([0 mm]);
title('y-acceptor & z-dark axis donor','fontsize',10);
% z-acceptor & z-dark axis donor
subplot(2,3,6)
pcolor(rhoa'*sin(phi),rhoa'*cos(phi),abs(xz).^2+abs(yz).^2); axis image; shading interp; set(gca,'fontsize',10)
xlabel('\itx\rm (nm)','fontsize',10)
ylabel('\ity\rm (nm)','fontsize',10)
caxis([0 mm]);
title('z-acceptor & z-dark axis donor','fontsize',10);

return

plot([-fliplr(rhoa) rhoa],[flipud(abs(yz(:,1)).^2+abs(zz(:,1)).^2); abs(yz(:,1)).^2+abs(zz(:,1)).^2]/max(abs(yz(:,1)).^2+abs(zz(:,1)).^2))
grid
print -dpng -r300 c:\Joerg\Doc\Microscopy\Defocused\Diamond\FRETMicroscopeAzDx
xlabel('\itx\rm (nm)'); ylabel('intensity (a.u.)'); title('y-acceptor & x-dark axis donor')
print -dpng -r300 c:\Joerg\Doc\Microscopy\Defocused\Diamond\FRETMicroscopeAzDx