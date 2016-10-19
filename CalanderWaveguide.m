% CalanderWaveguide

load metals
ng = 1.52;
nw = 1.333;
lam = 570;
nag = silver(wavelength==lam);
theta = (pi/2e4:pi/1e4:pi/2)';
dtheta = mean(diff(theta));
zv = 5:lam;
dv = 5:100;

nl = 2.5; % ref. index of spacer layer
clear vd1 pcd1 psd1 vu1 pcu1 psu1 lvd1 lvu1 lpd1 lpu1 qvd1 qvu1 qpd1 qpu1 vert1 vertup1 para1 paraup1
for j=1:length(zv)
    z = zv(j);
    for k= 1:length(dv)
        d = dv(k);
        [vd1,pcd1,psd1] = DipoleL(theta,0,[ng nag nl],nw,nw,[d z]/lam*2*pi,0,[]);
        [vu1,pcu1,psu1] = DipoleL(theta,0,nw,nw,[nl nag ng],[],0,[z d]/lam*2*pi);    
        [lvd1(j,k),lvu1(j,k),lpd1(j,k),lpu1(j,k),qvd1(j,k),qvu1(j,k),qpd1(j,k),qpu1(j,k)] = LifetimeL(0,[ng nag nl],nw,nw,[d z]/lam*2*pi,0,[]);
        vert1(j,k) = dtheta*sum(sin(theta).*abs(vd1).^2); %./(qvu1(j,k)+qvd1(j,k));
        vertup1(j,k) = dtheta*sum(sin(theta).*abs(vu1).^2); %./(qvu1(j,k)+qvd1(j,k));
        para1(j,k) = dtheta*sum(sin(theta).*(abs(pcd1).^2+abs(psd1).^2)); %./(qpu1(j,k)+qpd1(j,k))/2;
        paraup1(j,k) = dtheta*sum(sin(theta).*(abs(pcu1).^2+abs(psu1).^2)); %./(qpu1(j,k)+qpd1(j,k))/2;
    end
    j
end
save CalanderWaveguideAg25 zv dv lam ng nl nw vert1 vertup1 para1 paraup1 lvd1 lvu1 lpd1 lpu1 qvd1 qvu1 qpd1 qpu1


nl = 2.0; % ref. index of spacer layer
clear vd1 pcd1 psd1 vu1 pcu1 psu1 lvd1 lvu1 lpd1 lpu1 qvd1 qvu1 qpd1 qpu1 vert1 vertup1 para1 paraup1
for j=1:length(zv)
    z = zv(j);
    for k= 1:length(dv)
        d = dv(k);
        [vd1,pcd1,psd1] = DipoleL(theta,0,[ng nag nl],nw,nw,[d z]/lam*2*pi,0,[]);
        [vu1,pcu1,psu1] = DipoleL(theta,0,nw,nw,[nl nag ng],[],0,[z d]/lam*2*pi);    
        [lvd1(j,k),lvu1(j,k),lpd1(j,k),lpu1(j,k),qvd1(j,k),qvu1(j,k),qpd1(j,k),qpu1(j,k)] = LifetimeL(0,[ng nag nl],nw,nw,[d z]/lam*2*pi,0,[]);
        vert1(j,k) = dtheta*sum(sin(theta).*abs(vd1).^2); %./(qvu1(j,k)+qvd1(j,k));
        vertup1(j,k) = dtheta*sum(sin(theta).*abs(vu1).^2); %./(qvu1(j,k)+qvd1(j,k));
        para1(j,k) = dtheta*sum(sin(theta).*(abs(pcd1).^2+abs(psd1).^2)); %./(qpu1(j,k)+qpd1(j,k))/2;
        paraup1(j,k) = dtheta*sum(sin(theta).*(abs(pcu1).^2+abs(psu1).^2)); %./(qpu1(j,k)+qpd1(j,k))/2;
    end
    j
end
save CalanderWaveguideAg20 zv dv lam ng nl nw vert1 vertup1 para1 paraup1 lvd1 lvu1 lpd1 lpu1 qvd1 qvu1 qpd1 qpu1



nag = gold(wavelength==lam);

nl = 2.5; % ref. index of spacer layer
clear vd1 pcd1 psd1 vu1 pcu1 psu1 lvd1 lvu1 lpd1 lpu1 qvd1 qvu1 qpd1 qpu1 vert1 vertup1 para1 paraup1
for j=1:length(zv)
    z = zv(j);
    for k= 1:length(dv)
        d = dv(k);
        [vd1,pcd1,psd1] = DipoleL(theta,0,[ng nag nl],nw,nw,[d z]/lam*2*pi,0,[]);
        [vu1,pcu1,psu1] = DipoleL(theta,0,nw,nw,[nl nag ng],[],0,[z d]/lam*2*pi);    
        [lvd1(j,k),lvu1(j,k),lpd1(j,k),lpu1(j,k),qvd1(j,k),qvu1(j,k),qpd1(j,k),qpu1(j,k)] = LifetimeL(0,[ng nag nl],nw,nw,[d z]/lam*2*pi,0,[]);
        vert1(j,k) = dtheta*sum(sin(theta).*abs(vd1).^2); %./(qvu1(j,k)+qvd1(j,k));
        vertup1(j,k) = dtheta*sum(sin(theta).*abs(vu1).^2); %./(qvu1(j,k)+qvd1(j,k));
        para1(j,k) = dtheta*sum(sin(theta).*(abs(pcd1).^2+abs(psd1).^2)); %./(qpu1(j,k)+qpd1(j,k))/2;
        paraup1(j,k) = dtheta*sum(sin(theta).*(abs(pcu1).^2+abs(psu1).^2)); %./(qpu1(j,k)+qpd1(j,k))/2;
    end
    j
end
save CalanderWaveguideAu25 zv dv lam ng nl nw vert1 vertup1 para1 paraup1 lvd1 lvu1 lpd1 lpu1 qvd1 qvu1 qpd1 qpu1


nl = 2.0; % ref. index of spacer layer
clear vd1 pcd1 psd1 vu1 pcu1 psu1 lvd1 lvu1 lpd1 lpu1 qvd1 qvu1 qpd1 qpu1 vert1 vertup1 para1 paraup1
for j=1:length(zv)
    z = zv(j);
    for k= 1:length(dv)
        d = dv(k);
        [vd1,pcd1,psd1] = DipoleL(theta,0,[ng nag nl],nw,nw,[d z]/lam*2*pi,0,[]);
        [vu1,pcu1,psu1] = DipoleL(theta,0,nw,nw,[nl nag ng],[],0,[z d]/lam*2*pi);    
        [lvd1(j,k),lvu1(j,k),lpd1(j,k),lpu1(j,k),qvd1(j,k),qvu1(j,k),qpd1(j,k),qpu1(j,k)] = LifetimeL(0,[ng nag nl],nw,nw,[d z]/lam*2*pi,0,[]);
        vert1(j,k) = dtheta*sum(sin(theta).*abs(vd1).^2); %./(qvu1(j,k)+qvd1(j,k));
        vertup1(j,k) = dtheta*sum(sin(theta).*abs(vu1).^2); %./(qvu1(j,k)+qvd1(j,k));
        para1(j,k) = dtheta*sum(sin(theta).*(abs(pcd1).^2+abs(psd1).^2)); %./(qpu1(j,k)+qpd1(j,k))/2;
        paraup1(j,k) = dtheta*sum(sin(theta).*(abs(pcu1).^2+abs(psu1).^2)); %./(qpu1(j,k)+qpd1(j,k))/2;
    end
    j
end
save CalanderWaveguideAu20 zv dv lam ng nl nw vert1 vertup1 para1 paraup1 lvd1 lvu1 lpd1 lpu1 qvd1 qvu1 qpd1 qpu1

return 

para1 = para1/2;
paraup1 = paraup1/2;

subplot(141)
pcolor(dv,zv,vert1./(qvu1+qvd1));
caxis([0 1])
shading flat
hold on; contour(dv,zv,vert1./(qvu1+qvd1),[0:0.1:1],'linecolor','y','linewidth',1); hold off
ylabel('dielectric layer thickness [nm]')
title('glass')
subplot(142)
pcolor(dv,zv,vertup1./(qvu1+qvd1));
caxis([0 1])
shading flat
hold on; contour(dv,zv,vertup1./(qvu1+qvd1),[0:0.1:1],'y','linewidth',1); hold off
xlabel('metal film thickness [nm]');
set(gca,'ytick',[])
title('water')
subplot(143)
pcolor(dv,zv,(qvd1-vert1)./(qvu1+qvd1));
caxis([0 1])
shading flat
hold on; contour(dv,zv,(qvd1-vert1)./(qvu1+qvd1),[0:0.1:1],'y','linewidth',1); hold off
set(gca,'ytick',[])
title('metal')
subplot(144)
balken(0:0.1:1,'v');

subplot(141)
contourf(dv,zv,vert1./(qvu1+qvd1),[0:0.1:1],'k');
caxis([-0.1 1.])
ylabel('dielectric layer thickness [nm]')
title('glass')
subplot(142)
contourf(dv,zv,vertup1./(qvu1+qvd1),[0:0.1:1],'k');
caxis([-0.1 1.])
xlabel('metal film thickness [nm]');
set(gca,'ytick',[])
title('water')
subplot(143)
contourf(dv,zv,(qvd1-vert1)./(qvu1+qvd1),[0:0.1:1],'k');
caxis([-0.1 1.])
set(gca,'ytick',[])
title('metal')
subplot(144)
balken(0:0.1:1,'v');
caxis([-0.1 1.])
set(gca,'fontsize',16)

figure

subplot(141)
pcolor(dv,zv,para1./(qpu1+qpd1));
caxis([0 1])
shading flat
ylabel('dielectric layer thickness [nm]')
title('glass')
subplot(142)
pcolor(dv,zv,paraup1./(qpu1+qpd1));
caxis([0 1])
shading flat
xlabel('metal film thickness [nm]');
set(gca,'ytick',[])
title('water')
subplot(143)
pcolor(dv,zv,(qpd1-para1)./(qpu1+qpd1));
caxis([0 1])
shading flat
set(gca,'ytick',[])
title('metal')
subplot(144)
balken(0:0.1:1,'v');

subplot(141)
contourf(dv,zv,para1./(qpu1+qpd1),[0:0.1:1],'k');
caxis([-0.1 1])
ylabel('dielectric layer thickness [nm]')
title('glass')
subplot(142)
contourf(dv,zv,paraup1./(qpu1+qpd1),0:0.1:1,'k');
caxis([-0.1 1])
xlabel('metal film thickness [nm]');
set(gca,'ytick',[])
title('water')
subplot(143)
contourf(dv,zv,(qpd1-para1)./(qpu1+qpd1),0:0.1:1,'k');
caxis([-0.1 1])
set(gca,'ytick',[])
title('metal')
subplot(144)
balken(0:0.1:1,'v');
caxis([-0.1 1])
set(gca,'fontsize',16)