% Calander

load metals
ng = 1.52;
nw = 1.333;
lambda=501:700;
for j=1:200 [d(j), phi(j), rp(j)] = SPLocate([ng gold(wavelength==lambda(j)) nw],lambda(j)); end

[ax,h1,h2] = plotyy(lambda,d,lambda,phi/pi*180);
set(get(ax(1),'ylabel'),'string','layer thickness [nm]')
set(get(ax(2),'ylabel'),'string','incidence angle [°]')
xlabel('wavelength [nm]')
print -dpng -r600 'SP%Wavelength' 

lam = 570;
nag = gold(wavelength==lam);
theta = asin(nw/ng):pi/1000:pi/2;
w1 = ng*cos(theta);
dv = 10:0.5:90;
clear rp
for j=1:length(dv)
    rp(:,j) = fresnel(w1,[ng nag nw],dv(j)/lam*2*pi).';
end
contourf(dv,theta/pi*180,abs(rp).^2,[0:0.05:1]); xlabel('layer thickness [nm]'); ylabel('incidence angle [°]');
xlabel('layer thickness [nm]'); ylabel('incidence angle [°]');
caxis([0 1]);
colorbar
colormap jet
print -dpng -r600 'SP%Thickness&Angle_Col'
colormap(1-gray)
brighten(0.7)
print -dpng -r600 'SP%Thickness&Angle'


theta = pi/2000:pi/1000:pi/2;
zv = 2.5:2.5:lam;

clear v1 ps1 pc1 lvd1 vlu1 lpd1 lpu1 qvd1 qvu1 qpd1 qpu1
for j=1:length(zv)
    [vd1(:,j),pcd1(:,j),psd1(:,j)] = DipoleL(theta,zv(j)/lam*2*pi,[ng nag],nw,nw,d(lambda==lam)/lam*2*pi,zv(j)/lam*2*pi,[]);
    [vu1(:,j),pcu1(:,j),psu1(:,j)] = DipoleL(theta,0,nw,nw,[nag ng],[],zv(j)/lam*2*pi,d(lambda==lam)/lam*2*pi);    
    [lvd1(j),lvu1(j),lpd1(j),lpu1(j),qvd1(j),qvu1(j),qpd1(j),qpu1(j)] = LifetimeL(zv(j)/lam*2*pi,[ng nag],nw,nw,d(lambda==lam)/lam*2*pi,zv(j)/lam*2*pi,[]);
    j
end

dtheta = mean(diff(theta));

vert1 = dtheta*(sin(theta)*abs(vd1).^2)./(qvu1+qvd1);
vertup1 = dtheta*(sin(theta)*abs(vu1).^2)./(qvu1+qvd1);
ind = theta>asin(nw/ng);
vert1saf = dtheta*(sin(theta(ind))*abs(vd1(ind,:)).^2)./(qvu1+qvd1);
tst = vert1==max(vert1);
plot(zv,vert1,zv,vertup1,zv,(qvd1-lvd1)./(qvu1+qvd1))
hold on
plot(zv(tst)*[1 1],[0 1],':','linewidth',1.2)    
hold off
axis([0 lam 0 1])
xlabel('distance [nm]');
ylabel('partition of energy')
print -dpng -r600 VertSPCE

tmp=abs([flipud(vu1(:,tst)); vd1(:,tst)]).^2/(qvu1(tst)+qvd1(tst))/dtheta;
mpolar([pi-fliplr(theta) theta],tmp,1.1*max(tmp),[1,4],[ng nw])
print -dpng -r600 VertSPCEADR    

para1 = dtheta*(sin(theta)*(abs(pcd1).^2+abs(psd1).^2))./(qpu1+qpd1)/2;
paraup1 = dtheta*(sin(theta)*(abs(pcu1).^2+abs(psu1).^2))./(qpu1+qpd1)/2;
para1saf = dtheta*(sin(theta(ind))*(abs(pcd1(ind,:)).^2+abs(psd1(ind,:)).^2))./(qpu1+qpd1)/2;
tst = para1==max(para1);
plot(zv,para1,zv,paraup1,zv,(qpd1-lpd1)./(qpu1+qpd1))
hold on
plot(zv(tst)*[1 1],[0 1],':','linewidth',1.2)    
hold off
axis([0 lam 0 1])    
xlabel('distance [nm]');
ylabel('partition of energy')
print -dpng -r600 ParaSPCE

tmp=abs([flipud(pcu1(:,tst)+psu1(:,tst)); pcd1(:,tst)+psd1(:,tst)]).^2/(qvu1(tst)+qvd1(tst))/dtheta/2;
mpolar([pi-fliplr(theta) theta],tmp,1.1*max(tmp),[1,4],[ng nw])
print -dpng -r600 ParaSPCEADR

QY = [0.1 0.5 1];
S0 = 4/3*nw;
vert = qvu1+qvd1;
para = qpu1+qpd1;
for j=1:length(QY)
    kappa = sqrt(QY(j)*(vert-para)./(QY(j)*para+(1-QY(j))*S0));
    r = (QY(j)*vert+(1-QY(j))*S0)./(QY(j)*para+(1-QY(j))*S0);
    SPCE(:,j) = (1./(vert-para).*(vert1.*vert-para1.*para - (vert1.*vert-r.*para1.*para)./kappa.*atan(kappa)))';
    SPCEsaf(:,j) = (1./(vert-para).*(vert1saf.*vert-para1saf.*para - (vert1saf.*vert-r.*para1saf.*para)./kappa.*atan(kappa)))';
    SPCEvert(:,j) = (QY(j)*vert1.*vert./((1-QY(j))*S0+QY(j)*vert))';
    SPCEsafvert(:,j) = (QY(j)*vert1saf.*vert./((1-QY(j))*S0+QY(j)*vert))';
    SPCEpara(:,j) = (QY(j)*para1.*para./((1-QY(j))*S0+QY(j)*para))';
    SPCEsafpara(:,j) = (QY(j)*para1saf.*para./((1-QY(j))*S0+QY(j)*para))';    
end


% pure glass: dipole transmission, reflectivity, absorption + angular dist. of emission, lifetime

clear v2 ps2 pc2 lvd2 vlu2 lpd2 lpu2 qvd2 qvu2 qpd2 qpu2
zv = [0 zv];
for j=1:length(zv)
    [vd2(:,j),pcd2(:,j),psd2(:,j)] = DipoleL(theta,zv(j)/lam*2*pi,ng,nw,nw,[],zv(j)/lam*2*pi,[]);
    [vu2(:,j),pcu2(:,j),psu2(:,j)] = DipoleL(theta,0,nw,nw,ng,[],zv(j)/lam*2*pi,[]);    
    [lvd2(j),lvu2(j),lpd2(j),lpu2(j),qvd2(j),qvu2(j),qpd2(j),qpu2(j)] = LifetimeL(zv(j)/lam*2*pi,ng,nw,nw,[],zv(j)/lam*2*pi,[]);
    j
end

dtheta = mean(diff(theta));

vert2 = dtheta*(sin(theta)*abs(vd2).^2)./(qvu2+qvd2);
vertup2 = dtheta*(sin(theta)*abs(vu2).^2)./(qvu2+qvd2);
ind = theta>asin(nw/ng);
vert2saf = dtheta*(sin(theta(ind))*abs(vd2(ind,:)).^2)./(qvu2+qvd2);
plot(zv,vert2,zv,vertup2,zv,vert2saf)
axis([0 lam 0 1])
xlabel('distance [nm]');
ylabel('partition of energy')
print -dpng -r600 VertGlass

tst = 1;
tmp=abs([flipud(vu2(:,tst)); vd2(:,tst)]).^2/(qvu2(tst)+qvd2(tst))/dtheta;
mpolar([pi-fliplr(theta) theta],tmp,1.1*max(tmp),[1,4],[ng nw])
print -dpng -r600 VertGlassADR

para2 = mean(diff(theta))*(sin(theta)*(abs(pcd2).^2+abs(psd2).^2))./(qpu2+qpd2)/2;
paraup2 = dtheta*(sin(theta)*(abs(pcu2).^2+abs(psu2).^2))./(qpu2+qpd2)/2;
para2saf = dtheta*(sin(theta(ind))*(abs(pcd2(ind,:)).^2+abs(psd2(ind,:)).^2))./(qpu2+qpd2)/2;
plot(zv,para2,zv,paraup2,zv,para2saf)
axis([0 lam 0 1])    
xlabel('distance [nm]');
ylabel('partition of energy')
print -dpng -r600 ParaGlass

tmp=abs([flipud(pcu2(:,tst)+psu2(:,tst)); pcd2(:,tst)+psd2(:,tst)]).^2/(qvu2(tst)+qvd2(tst))/dtheta/2;
mpolar([pi-fliplr(theta) theta],tmp,1.1*max(tmp),[1,4],[ng nw])
print -dpng -r600 ParaGlassADR

plot(zv(2:end),4/3*nw./(qvu1+qvd1),':r',zv(2:end),4/3*nw./(qpu1+qpd1),':b',zv,4/3*nw./(qvu2+qvd2),'r',zv,4/3*nw./(qpu2+qpd2),'b')
ax = axis;
axis([0 lam 0 ax(4)])    
xlabel('distance [nm]');
ylabel('lifetime [\tau_0]')
print -dpng -r600 Lifetime

vert = qvu2+qvd2;
para = qpu2+qpd2;
for j=1:length(QY)
    kappa = sqrt(QY(j)*(vert-para)./(QY(j)*para+(1-QY(j))*S0));
    r = (QY(j)*vert+(1-QY(j))*S0)./(QY(j)*para+(1-QY(j))*S0);
    glass(:,j) = (1./(vert-para).*(vert2.*vert-para2.*para - (vert2.*vert-r.*para2.*para)./kappa.*atan(kappa)))';
    saf(:,j) = (1./(vert-para).*(vert2saf.*vert-para2saf.*para - (vert2saf.*vert-r.*para2saf.*para)./kappa.*atan(kappa)))';
    glassvert(:,j) = (QY(j)*vert2.*vert./((1-QY(j))*S0+QY(j)*vert))';
    safvert(:,j) = (QY(j)*vert2saf.*vert./((1-QY(j))*S0+QY(j)*vert))';
    glasspara(:,j) = (QY(j)*para2.*para./((1-QY(j))*S0+QY(j)*para))';
    safpara(:,j) = (QY(j)*para2saf.*para./((1-QY(j))*S0+QY(j)*para))';    
end

plot(zv,saf(:,1:2),zv,saf(:,3),'g',zv(2:end),SPCEsaf(:,1),':r',zv(2:end),SPCEsaf(:,2),':b',zv(2:end),SPCEsaf(:,3),':g')
ax = axis;
axis([0 lam 0 ax(4)])    
xlabel('distance [nm]');
ylabel('part of energy collected above critical angle')
print -dpng -r600 SAFCollectionEfficiency

plot(zv,glass(:,1:2),zv,glass(:,3),'g',zv(2:end),SPCE(:,1),':r',zv(2:end),SPCE(:,2),':b',zv(2:end),SPCE(:,3),':g')
ax = axis;
axis([0 lam 0 ax(4)])    
xlabel('distance [nm]');
ylabel('part of energy collected')    
print -dpng -r600 HalfspaceCollectionEfficiency

plot(zv,safvert(:,1:2),zv,safvert(:,3),'g',zv(2:end),SPCEsafvert(:,1),':r',zv(2:end),SPCEsafvert(:,2),':b',zv(2:end),SPCEsafvert(:,3),':g')
ax = axis;
axis([0 lam 0 ax(4)])    
xlabel('distance [nm]');
ylabel('part of energy collected above critical angle')
print -dpng -r600 SAFCollectionEfficiencyVert

plot(zv,safpara(:,1:2),zv,safpara(:,3),'g',zv(2:end),SPCEsafpara(:,1),':r',zv(2:end),SPCEsafpara(:,2),':b',zv(2:end),SPCEsafpara(:,3),':g')
ax = axis;
axis([0 lam 0 ax(4)])    
xlabel('distance [nm]');
ylabel('part of energy collected above critical angle')
print -dpng -r600 SAFCollectionEfficiencyPara

plot(zv,glassvert(:,1:2),zv,glassvert(:,3),'g',zv(2:end),SPCEvert(:,1),':r',zv(2:end),SPCEvert(:,2),':b',zv(2:end),SPCEvert(:,3),':g')
ax = axis;
axis([0 lam 0 ax(4)])    
xlabel('distance [nm]');
ylabel('part of energy collected')    
print -dpng -r600 HalfspaceCollectionEfficiencyVert

plot(zv,glasspara(:,1:2),zv,glasspara(:,3),'g',zv(2:end),SPCEpara(:,1),':r',zv(2:end),SPCEpara(:,2),':b',zv(2:end),SPCEpara(:,3),':g')
ax = axis;
axis([0 lam 0 ax(4)])    
xlabel('distance [nm]');
ylabel('part of energy collected')    
print -dpng -r600 HalfspaceCollectionEfficiencyPara

% photostability

S0 = 4/3*nw;
plot(zv,vert2.*(qvu2+qvd2)/S0,zv,para2.*(qpu2+qpd2)/S0,zv(2:end),vert1.*(qvu1+qvd1)/S0,':r',zv(2:end),para1.*(qpu1+qpd1)/S0,':b')
ax = axis;
axis([0 lam 0 ax(4)])    
xlabel('distance [nm]');
ylabel('detectable photon yield until photobleaching')
print -dpng -r600 Photostabilty    

plot(zv,vert2saf.*(qvu2+qvd2)/S0,zv,para2saf.*(qpu2+qpd2)/S0,zv(2:end),vert1saf.*(qvu1+qvd1)/S0,':r',zv(2:end),para1saf.*(qpu1+qpd1)/S0,':b')
ax = axis;
axis([0 lam 0 ax(4)])    
xlabel('distance [nm]');
ylabel('SAF photon yield until photobleaching')
print -dpng -r600 PhotostabiltySAF


% excitation

lamex = 532;
zv = 0:2.5:lam;
nag = gold(wavelength==lamex);
rp = fresnel(ng*cos(theta),[ng nag nw],d(lambda==lam)/lam*2*pi);
psi = theta(abs(rp)==min(abs(rp)));
rp = fresnel(nw,[nw nag ng],d(lambda==lam)/lam*2*pi);
top = abs(exp(-i*nw*2*pi/lam*zv) + rp*exp(i*nw*2*pi/lam*zv)).^2;
[bla,bla,tp] = fresnel(ng*cos(psi),[ng nag nw],d(lambda==lam)/lam*2*pi);
plasmon = abs(tp*exp(i*sqrt(nw^2-ng^2*sin(psi)^2)*2*pi/lam*zv)).^2;
[bla,bla,tp] = fresnel(ng,ng,nw);
wide = abs(tp*exp(i*nw*2*pi/lam*zv)).^2;    
[bla,bla,tp] = fresnel(sqrt(ng^2-nw^2),ng,nw);    
tirf = abs(tp*exp(i*0*2*pi/lam*zv)).^2;

plot(zv,wide',zv,tirf',zv(2:end),top(2:end)',zv(2:end),plasmon(2:end)')
% legend({'wide','tirf','top','plasmon'})    
ax = axis;
axis([0 lam 0 ax(4)])    
xlabel('distance [nm]');
ylabel('excitation intensity [a.u.]')
print -dpng -r600 Excitation

plot(zv,saf(:,3).*wide',zv,saf(:,3).*tirf',zv(2:end),SPCEsaf(:,3).*top(2:end)',zv(2:end),SPCEsaf(:,3).*plasmon(2:end)')
% legend({'wide','tirf','top','plasmon'})    
ax = axis;
axis([0 lam 0 ax(4)])    
xlabel('distance [nm]');
ylabel('SAF brightness [a.u.]')
print -dpng -r600 BrightnessSAF

plot(zv,glass(:,3).*wide',zv,glass(:,3).*tirf',zv(2:end),SPCE(:,3).*top(2:end)',zv(2:end),SPCE(:,3).*plasmon(2:end)')
% legend({'wide','tirf','top','plasmon'})    
ax = axis;
axis([0 lam 0 ax(4)])    
xlabel('distance [nm]');
ylabel('brightness [a.u.]')
print -dpng -r600 Brightness

