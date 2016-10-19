function [p1 p2] = Brunstein(fname)

close all

im = double(imread(fname));
[fx, fy] = gradient(im);
[cl,len] = Cluster(abs(fy)+abs(fx)>mConv2(abs(fy)+abs(fx),disk(20)));
[x,y] = meshgrid(1:size(im,2),1:size(im,1));

tmp = cl==length(len);
x1 = sum(sum(x.*tmp))/sum(tmp(:));
y1 = sum(sum(y.*tmp))/sum(tmp(:));
[h,bin] = mHist(sqrt((x(:)-x1).^2+(y(:)-y1).^2),[],tmp(:));
r1 = bin*h/sum(h);
v1 = sqrt((bin-r1).^2*h/sum(h));
p1 = Simplex('RingFit',[x1 y1 r1 v1],[],[],[],[],tmp);
p1 = Simplex('RingFit',p1,[0 0 0 1],[inf inf inf 1],[],[],tmp);
p1 = Simplex('RingFit',p1,[0 0 0 1],[inf inf inf 1],[],[],tmp);

tmp = cl==length(len)-1;
x2 = sum(sum(x.*tmp))/sum(tmp(:));
y2 = sum(sum(y.*tmp))/sum(tmp(:));
[h,bin] = mHist(sqrt((x(:)-x2).^2+(y(:)-y2).^2),[],tmp(:));
r2 = bin*h/sum(h);
v2 = sqrt((bin-r2).^2*h/sum(h));
p2 = Simplex('RingFit',[x2 y2 r2 v2],[],[],[],[],tmp);
p2 = Simplex('RingFit',p2,[0 0 0 1],[inf inf inf 1],[],[],tmp);
p2 = Simplex('RingFit',p2,[0 0 0 1],[inf inf inf 1],[],[],tmp);

phi = 0:0.01:2*pi;
mim(im); 
hold on; 
plot(p1(1)+(p1(3)+p1(4)/2)*cos(phi),p1(2)+(p1(3)+p1(4)/2)*sin(phi),'c'); 
plot(p2(1)+(p2(3)+p2(4)/2)*cos(phi),p2(2)+(p2(3)+p2(4)/2)*sin(phi),'c'); 
hold off
text(3/4*size(im,2),20,['N.A. = ' mnum2str(p1(3)/p2(3),1,2)],'color','y')