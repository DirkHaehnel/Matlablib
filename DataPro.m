function [t1, t2, bs, bh, bsd, trans] = datapro(x);

wndw = 15;
offset = 150;
uniform = 1;
m = 2;
sig0 = 10;
mm = 10; 
fac = 3;

close all
x = x(:)';
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

bs = zeros(1, length(t1));
bh = zeros(1, length(t1));
for k = 1:length(t1)
    bs(k) = sum(y(t1(k):t2(k))-bg(t1(k):t2(k)));
    bh(k) = max(x(t1(k):t2(k)));
end
trans = (y(t1-1)-fac*bg(t1-1))./(fac*(bg(t1)-bg(t1-1))-y(t1)+y(t1-1)) + (t2-t1) + (y(t2)-fac*bg(t2))./(fac*(bg(t2+1)-bg(t2))-y(t2+1)+y(t2));

bs(t1<wndw-1) = [];
bh(t1<wndw-1) = [];
t2(t1<wndw-1) = [];
trans(t1<wndw-1) = [];
t1(t1<wndw-1) = [];
bs(t2>length(x)-wndw+1) = [];
bh(t2>length(x)-wndw+1) = [];
t1(t2>length(x)-wndw+1) = [];
trans(t2>length(x)-wndw+1) = [];
t2(t2>length(x)-wndw+1) = [];

bsd = hist(bs,1:max(bs));
bar(bsd); xlabel('Burst Size'); ylabel('Frequency'); figure
[tmp, ind] = sort(bs);
z = [];
j0=0;
for j=j0+1:j0+10
   tmp = round(0.5*(t1(ind(end-j+1))+t2(ind(end-j+1))));
   z = [z; x(tmp-wndw+1:tmp+wndw-1)];
end
plot(z'+offset*(cumsum(ones(size(z)))'-1),'k'); xlabel('Time'); ylabel('Photon Counts')
hold on
for j=j0+1:j0+10
    tmp = round(0.5*(t1(ind(end-j+1))+t2(ind(end-j+1))));
	plot((t1(ind(end-j+1)):t2(ind(end-j+1)))-tmp+wndw, x(t1(ind(end-j+1)):t2(ind(end-j+1))) + (j-j0-1)*offset,'ok');
end
hold off
figure
plot(bs,trans,'o','markersize',3); xlabel('Burst Size'); ylabel('Transit Time');

% Useful stuff
% for j=1:100 plot(t1(j)-5:t2(j)+5,x(t1(j)-5:t2(j)+5),t1(j):t2(j),x(t1(j):t2(j)),'o',t1(j)-5:t2(j)+5,bg(t1(j)-5:t2(j)+5)); pause; end
% m = mean(x); t=0:max(x); z=hist(x,t); zz=exp(t*log(m)-m-log(cumprod([1 1:max(t)])));
% tst=autocorr(x); semilogx((1:10000)*3e-3,tst(2:10001)); xlabel('Time (sec)'); ylabel('Autocorrelation')