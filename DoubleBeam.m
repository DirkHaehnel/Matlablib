% program for calculating the absolute brightness of single molecules for double beam excitation

wndw = 0;
offset = 150;
uniform = 1;
m = 2;
sig0 = 10;
mm = 10; 
fac = 2;

if exist('x1')==0 | exist('x2')==0
	disp('There is no data containing variables x1, x2'); break;
end
close all
x1 = x1(:)';
x2 = x2(:)';
x = x1 + x2;
mn = conv(x, ones(1,2*m+1)); 
mn = mn(m+1:end-m);
mn = mn./[(m+1):(2*m) (2*m+1)*ones(1,length(mn)-2*m) (2*m):-1:(m+1)]; 
sig = conv((x - mn).^2, ones(1,2*m+1)); 
sig = sig(m+1:end-m);
sig = sig./[(m+1):(2*m) (2*m+1)*ones(1,length(sig)-2*m) (2*m):-1:(m+1)]; 
y = mn + (x - mn).*sig./(sig + sig0);

if uniform
	bg = mean(x)*ones(size(x));
else
	bg = conv(x, ones(1,2*mm+1));
	bg = bg(mm+1:end-mm);
	bg = bg./[(mm+1):(2*mm) (2*mm+1)*ones(1,length(bg)-2*mm) (2*mm):-1:(mm+1)]; 
end

t = 1:length(y);
t1=t(y(1:end-1)<=fac*bg(1:end-1) & y(2:end)>fac*bg(2:end)) + 1;
t2=t(y(2:end)<=fac*bg(2:end) & y(1:end-1)>fac*bg(1:end-1));
if length(t2)==length(t1)+1 
   t2(1)=[]; 
end
if length(t1)==length(t2)+1 
   t1(end)=[]; 
end

bs1 = zeros(1, length(t1));
bs2 = bs1;
for k = 1:length(t1)
   bs1(k) = sum(x1(t1(k):t2(k)));
   bs2(k) = sum(x2(t1(k):t2(k)));   
end

bx = [t(1) reshape([t1;t1;t2;t2],1,4*length(t1)) t(end)];
by = [0 reshape(repmat([0;1;1;0],length(t1),1),1,4*length(t1)) 0];
plot(bx,by+48,t,x1)
