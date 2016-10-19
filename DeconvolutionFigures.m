% DeconvolutionFigures

clear all
close all

nref = 1.33;
NA = 1.2;
smax = NA/nref;

[q,w] = meshgrid(0.0025:0.005:4,-4:0.005:4);
amp0 = zeros(size(q));
amp0(q.^2+w.^2>0.98 & q.^2+w.^2<1.02 & w<0 & abs(atan(q./w))<asin(1.2/1.33))=1;

qq = 0:0.5:2e3; nmax = numel(qq);
qv = q(1,:)';
btrans = zeros(numel(qv),nmax);
for j=1:nmax
    btrans(:,j) = besselj(0,qv/max(qv)*qq(j));
end

ff = amp0*(btrans.*repmat(qv,[1 nmax]));
ff = mConv(ff,flipud(ff));
amp1 = ff*(btrans.*repmat(qq,[numel(qv) 1]))';

amp1 = amp1./max(amp1(:));
amp1(amp1<1e-6) = 0;

ff = amp1*(btrans.*repmat(qv,[1 nmax]));
ff = mConv(ff,ff);
amp2 = ff*(btrans.*repmat(qq,[numel(qv) 1]))';

amp2 = amp2./max(amp2(:));
amp2(amp2<1e-6) = 0;

CombineImages(log10(cat(3,[fliplr(amp1) amp1],[fliplr(amp1) amp1],[fliplr(amp2) amp2])),1,3,[],{{'(A)',''} ,{'(B)',''},{'(C)',''}})
h = get(gca,'children');
for j=1:3
    set(h(j),'fontname','calibri','fontsize',8)
end
line(2*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
line(4*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
hold on
plot(2*size(q,2),size(q,1)*0.5,'xy','markersize',15,'linewidth',1);
plot(2*size(q,2),size(q,1)*0.5,'oy','markersize',15,'linewidth',1);
hold off
text(4*size(q,2),0.47*size(q,1),'=','color','yellow','fontname','symbol','fontsize',24,'horizontalalignment','center')
line(size(q,2)*[1 1]*0.2, [0.9 0.8]*size(q,1),'color','yellow','linewidth',1)
line([0.2 0.4]*size(q,2),0.9*size(q,1)*[1 1],'color','yellow','linewidth',1)
text(0.46*size(q,2),0.9*size(q,1),'k_x','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic')
text(0.2*size(q,2),0.73*size(q,1),'k_z','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic','HorizontalAlignment','center')
caxis([-4 0])
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',7)
set(h,'PlotBoxAspectRatio',[1 10 1])

amp3 = amp1>0;
ff = amp3*(btrans.*repmat(qv,[1 nmax]));
ff = mConv(ff,ff);
amp3 = ff*(btrans.*repmat(qq,[numel(qv) 1]))';

amp3 = amp3./max(amp3(:));
amp3(amp3<1e-6) = 0;

figure
CombineImages(log10(cat(3,[fliplr(amp1) amp1]>0,[fliplr(amp1) amp1]>0,[fliplr(amp3) amp3])),1,3,[],{{'(A)',''} ,{'(B)',''},{'(C)',''}})
h = get(gca,'children');
for j=1:3
    set(h(j),'fontname','calibri','fontsize',8)
end
line(2*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
line(4*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
hold on
plot(2*size(q,2),size(q,1)*0.5,'xy','markersize',15,'linewidth',1);
plot(2*size(q,2),size(q,1)*0.5,'oy','markersize',15,'linewidth',1);
hold off
text(4*size(q,2),0.47*size(q,1),'=','color','yellow','fontname','symbol','fontsize',24,'horizontalalignment','center')
line(size(q,2)*[1 1]*0.2, [0.9 0.8]*size(q,1),'color','yellow','linewidth',1)
line([0.2 0.4]*size(q,2),0.9*size(q,1)*[1 1],'color','yellow','linewidth',1)
text(0.46*size(q,2),0.9*size(q,1),'k_x','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic')
text(0.2*size(q,2),0.73*size(q,1),'k_z','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic','HorizontalAlignment','center')
caxis([-4 0])
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',7)
set(h,'PlotBoxAspectRatio',[1 10 1])

print -r600 -dpng DeconvFigOTF

[x,y] = meshgrid(0.0025:0.005:2,-3:0.005:3);
xv = x(1,:)*2*pi;
yv = y(:,1)*2*pi;

int0 = exp(1i*yv*w(:,1)')*amp0*(repmat(qv,[1 numel(xv)]).*besselj(0,qv*xv));
int0 = int0/max(abs(int0(:)));
figure; 
mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,abs([fliplr(int0) int0].^2)); 
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.145 1])


int1 = real(exp(1i*yv*w(:,1)')*amp1*(repmat(qv,[1 numel(xv)]).*besselj(0,qv*xv)));
int1 = int1/max(abs(int1(:)));
figure; 
mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,real([fliplr(int1) int1])); 
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.145 1])


figure; 
mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,abs([fliplr(int0) int0].^4)); 
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.145 1])


int2 = real(exp(1i*yv*w(:,1)')*amp2*(repmat(qv,[1 numel(xv)]).*besselj(0,qv*xv)));
int2 = int2/max(abs(int2(:)));
figure; 
mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,real([fliplr(int2) int2])); 
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.145 1])


int3 = real(exp(1i*yv*w(:,1)')*amp3*(repmat(qv,[1 numel(xv)]).*besselj(0,qv*xv)));
int3 = int3/max(abs(int3(:)));
figure; 
mpcolor([-fliplr(xv) xv]/2/pi,yv'/2/pi,real([fliplr(int3) int3])); 
set(gca,'xtick',[],'ytick',[]);
line([-1.6 -0.6], [-2.6 -2.6], 'color', 'yellow', 'linewidth',2)
line([-1.6 -1.6], [-2.6 -1.6], 'color', 'yellow', 'linewidth',2)
text(-0.425,-2.55,'x','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
text(-1.6,-1.3,'z','color','yellow','fontname','calibri','fontsize',20,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',20,'OuterPosition',[0.75 0 0.145 1])

figure
dx = mean(diff(xv));
dy = mean(diff(yv));
CombineImages(cat(3,[fliplr(int1) int1],[fliplr(int2) int2],[fliplr(int3) int3]),1,3,[],{{'(A)',''} ,{'(B)',''},{'(C)',''}})
h = get(gca,'children');
for j=1:3
    set(h(j),'fontname','calibri','fontsize',12)
end
line([80 80+2*pi/dx], (size(int1,1)-80)*[1 1], 'color', 'yellow', 'linewidth',1)
line(80*[1 1], (size(int1,1)-80)+[0 -2*pi/dy], 'color', 'yellow', 'linewidth',1)
text(130+2*pi/dx,size(int1,1)-88,'x','color','yellow','fontname','calibri','fontsize',12,'fontangle','italic','HorizontalAlignment','center')
text(85,size(int1,1)-145-2*pi/dy,'z','color','yellow','fontname','calibri','fontsize',12,'fontangle','italic','HorizontalAlignment','center')
h = colorbar;
set(h,'LineWidth',0.5)
set(h,'Fontname','calibri','FontSize',11)

print -dpng -r600 DeconvFigPSF

