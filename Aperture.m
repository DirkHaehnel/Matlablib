function NA = Aperture(im)

NAmax = 1.65; 

[fx,fy] = gradient(im);

fx = fx./(abs(im)+(im==0));
fy = fy./(abs(im)+(im==0));

x = 1:size(im,2);
b1 = x(sum(fx(round(end/2-end/4:end/2+end/4),:))==max(sum(fx(round(end/2-end/4:end/2+end/4),1:round(end/4)))));
b1 = min(b1);
b2 = x(sum(fx(round(end/2-end/4:end/2+end/4),:))==min(sum(fx(round(end/2-end/4:end/2+end/4),round(3*end/4):end))));
b2 = max(b2);
x = 1:size(im,1);
b3 = x(sum(fy(:,round(end/2-end/4:end/2+end/4))')==max(sum(fy(1:round(end/4),round(end/2-end/4:end/2+end/4))')));
b3 = min(b3);
b4 = x(sum(fy(:,round(end/2-end/4:end/2+end/4))')==min(sum(fy(round(3*end/4):end,round(end/2-end/4:end/2+end/4))')));
b4 = max( b4);

rad1 = (b2-b1+b4-b3)/4;

x  = 3+(b4-b3)/2;
y  = 3+(b2-b1)/2;

e = im(b3-2:b4+2,b1-2:b2+2);

[yy,xx] = meshgrid(1:size(e,2),1:size(e,1));

rr = sqrt((xx-x).^2 + (yy-y).^2);

rd = 1:max(rr(:));
tmp = mhist(rr(:),rd);
s = mhist(rr(:),rd,e(:))./(tmp+(tmp==0));
s = s/max(s);

zv = linspace(0,0.1*pi,100);
n0 = 1.51;
n = 1.51;
n1 = 1.33;
d0 = [];
d = max(zv);
d1 = [];
theta = (0:pi/2e3:pi/2-1e-10)';
cc = 1./sqrt(cos(theta));
rm = n0*sin(theta);
intx = zeros(numel(theta),1);
inty = intx;
intz = intx;
for j=1:length(zv)
    [v, pc, ps] = DipoleL(theta,zv(j),n0,n,n1,d0,d,d1);
    intx = intx + cc.*abs(pc).^2;
    inty = inty + cc.*abs(ps).^2;
    intz = intz + cc.*abs(v).^2;
end

rdmax = find(s==max(s));

p = simplex('LuruFit1',[rdmax 1],[0 1],[inf 10],[],[],rd(rdmax-40:rdmax+20),s(rdmax-40:rdmax+20),rm,[intx+inty intz]);

rad2 = p(1);

mim(e);
hold on
phi = 0:pi/180:2*pi;
plot(y,x,'xw')
plot(y+rad1*cos(phi),x+rad1*sin(phi),'w');
plot(y+rad2*cos(phi),x+rad2*sin(phi),'c');
hold off
drawnow

NA = rad1/rad2
