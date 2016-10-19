% program HongPro for analyzing Hong's exonuclease cleavage data

global bld dt th

name = 'hh090297kk';
eval(['load ', name, '.txt']);
eval(['x = ', name, ';']);
eval(['clear ', name]);
close all

tt = 5:200;
t = x(tt,1);
% y = x(tt,[2 4:6]); c = [20 10 5 2.5]; n = 26; th = 0.01; p = [0.0084 0.1045 0.5821];
% y = x(tt,[13 15:17]); c = [20 10 5 2.5]; n = 26; th = 0.01; p = [0.0084 0.1045 0.5821];
% y = x(tt,[18:19 21:23]); c = [40 20 10 5 2.5]; n = 16; th = 0; p = [0.0084 0.1045 0.5821]; 
y = x(tt,[18:19 21]); c = [40 20 10]; n = 16; th = 0; p = [0.0084 0.1045 0.5821]; 
% y = x(tt,[7:9 11:12]); c = [40 20 10 5 2.5]; n = 16; th = 0.01; p = [0.0084 0.1045 0.5821]; 

dt = 0.1

for k=1:size(y,2)
   tmp=conv(y(:,k),[-1 2 -1]); tmp(1)=[]; tmp(end)=[]; tmp = tmp>0.3; 
   ind = 1:size(y,1); ind = ind(tmp); ind(ind==1)=[]; ind(ind==size(y,1))=[];
   y(ind,k) = 0.5*(y(ind-1,k)+y(ind+1,k));
end

close;
hold on
plot(t, y, 'o', 'MarkerSize', 2, 'EraseMode', 'none');
for k=1:length(c)
   bld(k) = plot(t, y(:,k), 'EraseMode','xor');
end
drawnow

p = SIMPLEX('hongfun', p, [1e-6 0 0], [], [], [], t, y, c, n);
   
kon = p(1);
koff = p(2);
kcl = p(3);
% coef = (exp((1-p(4))/p(5))+1)./(exp(((1:n)-p(4))/p(5))+1);
coef = ones(1,n);
f = zeros(size(y));
for j=1:length(c)
	z = zeros(2*n, 1);
	z(1) = 1;
	M = [(1-dt*c(j)*kon)*eye(n) dt*koff*eye(n); dt*c(j)*kon*eye(n) (1-dt*(koff+kcl))*eye(n) + dt*kcl*diag(ones(1,n-1),-1)];
   z = (M^round(t(1)/dt))*z;
   f(1,j) = sum(z);
	M = M^round((t(2)-t(1))/dt);
	for k=2:length(t)
      z = M*z;
      f(k,j) = [coef coef]*z;
   end;
   ind = y(:,j)>th;
   f(:,j) = [f(:,j) ones(length(t),1)]*abs([f(ind,j) ones(sum(ind),1)]\y(ind,j));
end
hold off
% plot(t,y(:,1),'o',t,y(:,2),'*',t,y(:,3),'s',t,y(:,4),'d',t,y(:,5),'v','MarkerSize',5);
% plot(t,y(:,1),'o',t,y(:,2),'*',t,y(:,3),'s',t,y(:,4),'d','MarkerSize',5);
plot(t,y(:,1),'o',t,y(:,2),'*',t,y(:,3),'s','MarkerSize',5);
axis([0 max(t) 0 1.4]);
hold on
plot(t,f,'LineWidth',1);
hold off
text(0.7*max(t), 1.3, ['k_{on} = ', sprintf('%0.3g',p(1))]);
text(0.7*max(t), 1.15, ['k_{off} = ', sprintf('%0.3g',p(2))]);
text(0.7*max(t), 1.0, ['k_{cl} = ', sprintf('%0.3g',p(3))]);
times16
xlabel('Time (sec)','FontName','Times'); ylabel('rel. Fluorescence','FontName','Times');
