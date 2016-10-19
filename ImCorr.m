function [y1,y2,c,p] = ImCorr(data1,data2,width,flag);

% [y1,y2,c,p] = ImCorr(data1,data2,width) correlates two images data1 and data2, where the optional parameter 
% width is the assumed maximum shift of one image relative to the other. Defaults value of width = 15.
% Output variables are the alligned image data y1 and y2; c gives the image correlation function. 
% The parameter p gives the x- and y-sigma values of the correlation distribution.
%
% Creator: Joerg Enderlein, 2001

if nargin<3 | isempty(width)
   width = 15;
end

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

if nargin==4
   imagesc(data1); drawnow
   [a,b] = ginput(2);
   a = round(a);
   b = round(b);   
   x1 = data1(b(1):b(2),a(1):a(2));
   x2 = data2(b(1):b(2),a(1):a(2));
   [n,m] = size(x1);
else
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
sj = round(sign(jm)/2+0.25);
sk = round(sign(km)/2+0.25);      
y1 = data1((1+sj*jm):(end+jm*(1-sj)),(1+sk*km):(end+km*(1-sk)));
y2 = data2((1-(1-sj)*jm):(end-jm*sj),(1-(1-sk)*km):(end-km*sk));
p = [sj sk];
