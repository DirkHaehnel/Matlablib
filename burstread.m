% Program burstread for finding bursts in SPC-data

zchar = name;
for j=1:floor(log10(num))+1 zchar = [zchar sprintf('%c', '0')]; end
[x, dtau] = spc430read(zchar);
bins = size(x,2);
channels = size(x,1);
x = sum(x);
for j=1:num-1
   zchar = num2str(j); 
   for j=1:floor(log10(num))-floor(log10(j)) zchar = [sprintf('%c', '0') zchar]; end
   disp(zchar);
   eval(['x = [x sum(spc430read(' char(39) name zchar char(39) '))];']);
end

m = 2;
mn = conv(x, ones(1,2*m+1)); 
mn = mn(m+1:end-m);
mn = mn./[(m+1):(2*m) (2*m+1)*ones(1,length(mn)-2*m) (2*m):-1:(m+1)]; 
sig = conv((x - mn).^2, ones(1,2*m+1)); 
sig = sig(m+1:end-m);
sig = sig./[(m+1):(2*m) (2*m+1)*ones(1,length(sig)-2*m) (2*m):-1:(m+1)]; 
sig0 = mean(sig);
y = mn + (x - mn).*sig./(sig + sig0);

bg = 3*mean(x);

t = 1:length(y);
t1=t(y(1:end-1)<=bg & y(2:end)>bg) + 1;
t2=t(y(2:end)<=bg & y(1:end-1)>bg);
if length(t2)==length(t1)+1 
   t2(1)=[]; 
end

burst = zeros(channels,length(t1));
zchar = name;
for j=1:floor(log10(num))+1 zchar = [zchar sprintf('%c', '0')]; end
x = spc430read(zchar);
zchar(end) = '1';
x = [x spc430read(zchar)];
cnt = 1;
for j=1:length(t1)
   disp(j);
   if t2(j)>=2*bins
      x(:,1:bins) = x(:,bins+1:end);
      cnt = cnt + 1;
	   zchar = num2str(cnt); 
	   for k=1:floor(log10(num))-floor(log10(cnt)) zchar = [sprintf('%c', '0') zchar]; end
      eval(['x(:,bins+1:end) = spc430read(' char(39) name zchar char(39) ');']);
      t1 = t1-bins;
      t2 = t2-bins;
   end
   burst(:,j) = sum(x(:,t1(j):t2(j))')';
end

burst(:,std(burst)==0)=[];
c=polyfit(5:40,log(sum(burst(5:40,:)')-min(sum(burst'))),1);
semilogy(1:channels,sum(burst')-min(sum(burst')),1:channels,sum(burst'),1:channels,exp(polyval(c,1:channels)));
drawnow;
tau=1/abs(c(1))*dtau;

eval(['save ' char(39) name char(39) ' burst dtau tau'])

