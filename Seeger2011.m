% program Seeger paper Opt Ex 2011
clear all
close all

theta = pi/1e4:pi/5e3:pi/2;
zv = 0:pi/20:2*pi;
n0 = 1.523;
n = 1.38;%1.333;
d0 = [];
d = 0;
d1 = [];
thetacut = [62 66 70]/180*pi;

for jz=1:length(zv)
    
[vd,pcd,psd] = DipoleL(theta,zv(jz),n0,n,n,d0,d,d1);
[vu,pcu,psu] = DipoleL(theta,1e-10,n,n,n0,d0,zv(jz),d1);

vtotal(jz) = pi/1e3*sin(theta)*(abs(vd).^2+abs(vu).^2);
ptotal(jz) = 0.5*pi/1e3*sin(theta)*(abs(pcd).^2+abs(pcu).^2+abs(psd).^2+abs(psu).^2);

for jt=1:length(thetacut)
    ind = theta>thetacut(jt) & theta<80/180*pi;
    vrel(jz,jt) = pi/1e3*sin(theta(ind))*abs(vd(ind)).^2;
    prel(jz,jt) = 0.5*pi/1e3*sin(theta(ind))*(abs(pcd(ind)).^2+abs(psd(ind)).^2);
end

jz

end

ptotal = ptotal'*ones(1,size(vrel,2));
vtotal = vtotal'*ones(1,size(vrel,2));
dtotal = sqrt(ptotal-vtotal);
res = abs((sqrt(ptotal).*(prel-vrel).*dtotal+(ptotal.*vrel-vtotal.*prel).*atanh(dtotal./sqrt(ptotal)))./...
    sqrt(ptotal)./dtotal.^3);

plot(zv/2/pi,res)