psi = linspace(0,asin(1.2/1.33),1e2);
[x,y] = meshgrid(0.0025:0.005:2,-3:0.005:3);
xv = x(1,:)*2*pi;
yv = y(:,1)*2*pi;
int0 = zeros(size(x));
for j=1:numel(psi)
    int0 = int0 + exp(1i*cos(psi(j))*yv)*besselj(0,sin(psi(j))*xv)*sin(psi(j));
end
int0 = int0/max(abs(int0(:)));

mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,[fliplr(real(int0)) real(int0)])
% xlabel('\rho (\lambda)','Fontname','Calibri','FontSize',20);
% ylabel('\itz\rm (\lambda)','Fontname','Calibri','FontSize',20)
%set(gca,'xticklabel',abs(get(gca,'xtick')),'Fontname','Calibri','FontSize',20)
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.145 1])

print -r600 -dpng PSFWidefieldA

mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,abs([fliplr(int0) int0]).^2)
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.13 1])

print -r600 -dpng PSFWidefieldB

mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,abs([fliplr(int0) int0]).^4)
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.13 1])

print -r600 -dpng PSFConfocal

mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,real([fliplr(int0) int0]).^2)
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.13 1])

print -r600 -dpng PSF4PiExc

mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,real([fliplr(int0) int0]).^2.*abs([fliplr(int0) int0]).^2)
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.13 1])

print -r600 -dpng PSF4PiA

mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,real([fliplr(int0) int0]).^4)
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.13 1])
set(h,'ytick',0.2:0.2:1)

print -r600 -dpng PSF4PiC

psi = linspace(asin(1.15/1.33),asin(1.2/1.33),1e1);
int1 = zeros(size(x));
for j=1:numel(psi)
    int1 = int1 + exp(1i*cos(psi(j))*yv)*besselj(0,sin(psi(j))*xv)*sin(psi(j));
end
int1 = int1/max(abs(int1(:)));

mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,abs([fliplr(int1) int1]).^2)
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.13 1])
set(h,'ytick',0.2:0.2:1)

print -r600 -dpng PSFSIMExc

mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,abs([fliplr(int1) int1]).^2.*abs([fliplr(int0) int0]).^2)
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.13 1])
set(h,'ytick',0.2:0.2:1)

print -r600 -dpng PSFSIM

