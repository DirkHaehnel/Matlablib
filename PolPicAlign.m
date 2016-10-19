function [ta1a, ta2a, ta1b, ta2b, ti1a, ti2a, ti1b, ti2b, c] = PolPicAlign(tag1a, tag2a, tag1b, tag2b, tim1a, tim2a, tim1b, tim2b);

% [ta1a, ta2a, ta1b, ta2b, ti1a, ti2a, ti1b, ti2b, c] = PolPicAlign(tag1a, tag2a, tag1b, tag2b, tim1a, tim2a, tim1b, tim2b) 
% correlates polarization images % tag1a, tag2a, tag1b, tag2b, where the number refers to the detecton 
% channel (emission polarization) and % the letters a, b refer to measurement (exctation polarization).
% Output variables are the aligned images; c gives the image correlation function. 
%
% Creator: Joerg Enderlein, 2001

width = 15;

data1 = tag1a + tag2a;
data2 = tag1b + tag2b;

[n1,m1] = size(data1);
[n2,m2] = size(data2);

if n1<n2
   data2 = data2(1:n1,:);
elseif n2<n1
   data1 = data1(1:n2,:);
end
if m1<m2
   data2 = data2(:,1:m1);
elseif m2<m1
   data1 = data1(:,1:m2);
end

x1 = data1;
x2 = data2;
[n,m] = size(x1);
if n>200;
   x1 = x1(1:200,:);
   x2 = x2(1:200,:);
   n = 200;
end
if m>200;
   x1 = x1(:,1:200);
   x2 = x2(:,1:200);
   m = 200;
end
   
h = waitbar(0,'Correlating images');
c = zeros(2*width+1, 2*width+1);
for j = -width:width
   for k = -width:width
      sj = round(sign(j)/2+0.25);
      sk = round(sign(k)/2+0.25);      
      c(width+1+j,width+1+k) = ...
         sum(sum(x1((1+sj*j):(end+j*(1-sj)),(1+sk*k):(end+k*(1-sk))).*x2((1-(1-sj)*j):(end-j*sj),(1-(1-sk)*k):(end-k*sk))))/...
         (n-abs(j))/(m-abs(k));
   end
   waitbar((j+width+1)/(2*width+1));
end
close(h);

[k,j] = meshgrid(-width:width,-width:width);
jm = j(c==max(max(c)));
km = k(c==max(max(c)));
[kk,jj] = meshgrid(max([km-3 -width]):0.1:min([km+3 width]),max([jm-3 -width]):0.1:min([jm+3 width]));
cc = interp2(k,j,c,kk,jj,'cubic');
jjm = jj(cc==max(max(cc)));
kkm = kk(cc==max(max(cc)));
tmp = polyfit(jj(kk==kkm),-log(cc(kk==kkm)),2);
p = sqrt(1/tmp(1)/2);
tmp = polyfit(kk(jj==jjm),-log(cc(jj==jjm)),2);
p = [p sqrt(1/tmp(1)/2)];
sj = round(sign(jm)/2+0.25);
sk = round(sign(km)/2+0.25);      
ta1a = tag1a((1+sj*jm):(end+jm*(1-sj)),(1+sk*km):(end+km*(1-sk)));
ta2a = tag2a((1+sj*jm):(end+jm*(1-sj)),(1+sk*km):(end+km*(1-sk)));
ta1b = tag1b((1-(1-sj)*jm):(end-jm*sj),(1-(1-sk)*km):(end-km*sk));
ta2b = tag2b((1-(1-sj)*jm):(end-jm*sj),(1-(1-sk)*km):(end-km*sk));
ti1a = tim1a((1+sj*jm):(end+jm*(1-sj)),(1+sk*km):(end+km*(1-sk)));
ti2a = tim2a((1+sj*jm):(end+jm*(1-sj)),(1+sk*km):(end+km*(1-sk)));
ti1b = tim1b((1-(1-sj)*jm):(end-jm*sj),(1-(1-sk)*km):(end-km*sk));
ti2b = tim2b((1-(1-sj)*jm):(end-jm*sj),(1-(1-sk)*km):(end-km*sk));

