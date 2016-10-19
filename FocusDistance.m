% FocusDistance

y = double(imread('D:\Joerg\Doc\Fcs\2Focus\ScannerCalibration\2006-08-28\Masstab.tif'));
close; p=simplex('GridFit',[70 100 10],[],[],[],[],y(:,100),2)
for j=2:4 close; p(:,j)=simplex('GridFit',p(:,j-1),[],[],[],[],y(:,j*100),2); end
c=polyfit(100:100:400,p(2,:),1);
pix=mean(10e3./(p(1,:)/cos(atan(c(1)))))
std(10e3./(p(1,:)/cos(atan(c(1)))))
% Maﬂstab = 132.95 ± 0.06 nm / pixel

x=double(imread('D:\Joerg\Doc\Fcs\2Focus\ScannerCalibration\2006-08-28\LinienLaser1linksLaser2rechts_2.tif'));
x = x(401:800,:);
a=1:size(x,1); t=1:size(x,2);
for j=1:size(x,2) r(j)=mean(a(x(:,j)==min(x(:,j)))); end
for j=1:length(r) close; rr(:,j)=simplex('Gauss',[r(j) 10],[],[],[],[],r(j)-20:r(j)+20,max(x(:))-x(r(j)-20:r(j)+20,j)); rr(:,j)=simplex('Gauss',rr(:,j),[],[],[],[],r(j)-20:r(j)+20,max(x(:))-x(r(j)-20:r(j)+20,j)); end
ind=t<400 | t>650;
tt=t(ind);
cc=[tt' double(tt<400)' double(tt>650)']\rr(1,ind)'
plot(t,rr(1,:),t,polyval(cc([1 2]),t),t,polyval(cc([1 3]),t))
(cc(3)-cc(2))*cos(atan(cc(1)))*pix*sqrt(2)

