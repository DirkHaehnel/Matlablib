% function gas()

close all
bdwidth = 5;
topwidth = 40;
set(0, 'units', 'pixels');
scnsize = get(0, 'screensize');
win = [bdwidth, bdwidth, scnsize(3)-2*bdwidth, scnsize(4)-(bdwidth+topwidth)];
figure('position', win);
N = 500;
col = ones(N,1);
row = col';
dt = 0.1;
p = N*rand(N,2);
v = 2*rand(N,1)-1; v = [v, sign(rand(N,1)-0.5).*sqrt(1-v.*v)];
subplot(121);
bld = plot(p(:,1), p(:,2), '.'); set(bld, 'erasemode', 'background'); axis image; 
axis([-1.5 N+1 -1 N+1.5]);
set(gca, 'xtick', [], 'ytick', []);
subplot(122);
t = 0.05:0.1:4; z = zeros(size(t));
hist(sqrt(sum((v.*v)')), t); axis([0 max(t) 0 N]);
axis square
drawnow;
for j=1:1e6
   p = p + dt*v;
   v(p(:,1)<0,1) = -v(p(:,1)<0,1);
   p(p(:,1)<0,1) = -p(p(:,1)<0,1);
   v(p(:,2)<0,2) = -v(p(:,2)<0,2);
   p(p(:,2)<0,2) = -p(p(:,2)<0,2);
   v(p(:,1)>N,1) = -v(p(:,1)>N,1);
   p(p(:,1)>N,1) = 2*N - p(p(:,1)>N,1);
   v(p(:,2)>N,2) = -v(p(:,2)>N,2);
   p(p(:,2)>N,2) = 2*N - p(p(:,2)>N,2);
   dist = (p(:,1)*row - col*p(:,1)').^2 + (p(:,2)*row - col*p(:,2)').^2;
   dist = dist>0 & dist<2*dt^2;
   
%    [ord, ord] = sort(p(:,1)+N*p(:,2));
%    p = p(ord,:); v = v(ord,:);
   d = diff(p(:,1))<dt & diff(p(:,2))<dt;
   if sum(d)>0
      dd = rand(sum(d),1); dd = [dd sqrt(1-dd.^2)];
      if sum(d)>1
         tmp = dd.*(sum(((v(logical([0; d]),:)-v(logical([d; 0]),:)).*dd)')'*ones(1,2)); 
      else
         tmp = dd.*((v(logical([0; d]),:)-v(logical([d; 0]),:)).*dd); 
      end         
      vold = v(logical([0; d]),:);
      v(logical([d; 0]),:) = v(logical([d; 0]),:) + tmp; 
      v(logical([0; d]),:) = vold - tmp;
   end
   if mod(j,100)==0 
      set(bld, 'xdata', p(:,1), 'ydata', p(:,2));
      r = hist(sqrt(sum((v.*v)')), t); bar(t,r); axis([0 max(t) 0 N/5]);
      % title(['average velocity = ' mnum2str(mean(sqrt(sum(v.*v,2))),2,3)]);
      axis square
      drawnow
   end
end

