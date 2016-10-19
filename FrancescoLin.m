function [tx, tz, hx, hy, hz, htx, hty, htz, gaussx, gaussy, gaussz, gausstx, gaussty, gausstz, kx, ky, kz, cx, cy, cz] = ...
   FrancescoLin(fname, flag);

% [tx, tz, hx, hy, hz, htx, hty, htz, gaussx, gaussy, gaussz, gausstx, gaussty, gausstz, kx, ky, kz, cx, cy, cz] = ...
% FrancescoLin(fname, flag)
% FrancescoLin evaluates Borwinian motion data of an optical tweezer
% fname - name of input data file
% if flag == 1, the results are saved on disk 
% t* - abscissa of histograms
% h* - histograms
% ht* - modeled Gaussians to histograms
% k* - tweezer stiffness
% c* - 4*diffusion/k*

close all
n = 100; % number of binnings
m = 10; % number of concatenated intervals
intake = 5e4; % how many data to process at once

if nargin==0
   [fname, pname] = uigetfile('*.*');
   fname = [pname fname];
   flag = questdlg('Save Results?');
   if ~(flag(1)=='y') flag=[]; end
   drawnow
end

fin = fopen(fname,'r');

an = max(union(0,union(findstr(fname,'\'),findstr(fname,'/'))))+1;
en = min(union(findstr(fname,'.'),length(fname)+1))-1;
fname = fname(an:en) ;

[data, count] = fscanf(fin, '%f ', intake);

% determine smalles voltage increment:
tmp = abs(diff(data)); tmp(tmp==0) = [];
dt = min(tmp);

% bulding up the histograms:
x = data(1:5:end) - data(2:5:end) - data(3:5:end) + data(4:5:end);
y = data(2:5:end) + data(1:5:end) - data(3:5:end) - data(4:5:end);
z = data(5:5:end);
w = abs(data(1:5:end) + data(2:5:end) + data(3:5:end) + data(4:5:end));

dx = x(n*m+1:end)-x(1:end-n*m);
dy = y(n*m+1:end)-y(1:end-n*m);
dz = z(n*m+1:end)-z(1:end-n*m);
tx = ceil(max([max(abs(dx)) max(abs(dy))])/dt);
tz = ceil(max(abs(dz))/dt);
tx = dt*(-tx:tx)/mean(w);
tz = dt*(-tz:tz);
x = x./w;
y = y./w;
hx = [];
hy = [];
hz = [];
for j = 1:n
   dx = x(j*m+1:end) - x(1:end-j*m);
   dy = y(j*m+1:end) - y(1:end-j*m);
   dz = z(j*m+1:end) - z(1:end-j*m);
   hx = [hx hist(dx, tx)'];
   hy = [hy hist(dy, tx)'];
   hz = [hz hist(dz, tz)'];
end
close all
scrsz = get(0,'ScreenSize'); 
figure('Name',['Histograms ' fname], 'NumberTitle','off', 'Position',[scrsz(3)/2-500 scrsz(4)/2-100 1000 400]);
subplot(131);
plot(tx,hx./(ones(size(hx,1),1)*sum(hx))); axis tight;
xlabel('X');
ylabel('probability');
subplot(132);
plot(tx,hy./(ones(size(hy,1),1)*sum(hy))); axis tight;
xlabel('Y');
subplot(133);

plot(tz,hz./(ones(size(hz,1),1)*sum(hz))); axis tight;
xlabel('Z');
drawnow
while count>0
   [data, count] = fscanf(fin, '%f ', intake);
   if count>0
      x = data(1:5:end) - data(2:5:end) - data(3:5:end) + data(4:5:end);
      y = data(2:5:end) + data(1:5:end) - data(3:5:end) - data(4:5:end);
      z = data(5:5:end);
      w = abs(data(1:5:end) + data(2:5:end) + data(3:5:end) + data(4:5:end));
      x = x./w;
      y = y./w;
      for j=1:n
         dx = x(j*m+1:end) - x(1:end-j*m);
         dy = y(j*m+1:end) - y(1:end-j*m);
         dz = z(j*m+1:end) - z(1:end-j*m);
         hx(:,j) = hx(:,j) + hist(dx, tx)';
         hy(:,j) = hy(:,j) + hist(dy, tx)';
         hz(:,j) = hz(:,j) + hist(dz, tz)';
      end
      subplot(131);
      plot(tx,hx./(ones(size(hx,1),1)*sum(hx))); axis tight;
      xlabel('X');
      ylabel('probability');
      subplot(132);
      plot(tx,hy./(ones(size(hy,1),1)*sum(hy))); axis tight;
      xlabel('Y');
      subplot(133);
      plot(tz,hz./(ones(size(hz,1),1)*sum(hz))); axis tight;
      xlabel('Z');
      drawnow
   end
end
fclose(fin);

% normalize distributions:
hx = hx./(ones(size(hx,1),1)*sum(hx));
hy = hy./(ones(size(hy,1),1)*sum(hy));
hz = hz./(ones(size(hz,1),1)*sum(hz));
hx = hx(2:end-1,:);
hy = hy(2:end-1,:);
hz = hz(2:end-1,:);
tx = tx(2:end-1);
tz = tz(2:end-1);

% plot final results
subplot(131);
plot(tx,hx);
xlabel('X');
ylabel('probability');
subplot(132);
plot(tx,hy);
xlabel('Y');
subplot(133);
plot(tz,hz);
xlabel('Z');
drawnow

% calculate Gaussian parameters and distributions:
ttx = tx'*ones(1,size(hx,2));
ttz = tz'*ones(1,size(hz,2));
mx = sum(ttx.*hx);
my = sum(ttx.*hy);
mz = sum(ttz.*hz);
dx = sqrt(sum((ttx-ones(size(ttx,1),1)*mx).^2.*hx));
dy = sqrt(sum((ttx-ones(size(ttx,1),1)*my).^2.*hy));
dz = sqrt(sum((ttz-ones(size(ttz,1),1)*mz).^2.*hz));
gaussx = []; gaussy = []; gaussz = [];
hh = figure('Name',['Gaussian fits ' fname], 'NumberTitle','off');
htx = []; hty = []; htz = [];
for j = 1:n
   tmp = simplex('GaussFun',[mx(j) dx(j)],[-inf 0],[],[],[],tx,hx(:,j));
   gaussx = [gaussx tmp(2)];
   [tmp, tmp] = GaussFun(tmp,tx,hx(:,j));
   htx = [htx tmp];
   tmp = simplex('GaussFun',[my(j) dy(j)],[-inf 0],[],[],[],tx,hy(:,j));
   gaussy = [gaussy tmp(2)];
   [tmp, tmp] = GaussFun(tmp,tx,hy(:,j));
   hty = [hty tmp];
   tmp = simplex('GaussFun',[mz(j) dz(j)],[-inf 0],[],[],[],tz,hz(:,j));
   gaussz = [gaussz tmp(2)];   
   [tmp, tmp] = GaussFun(tmp,tz,hz(:,j));
   htz = [htz tmp];
end
close(hh);

% calculate final parameters:
figure('Name',['Parameter fits ' fname], 'NumberTitle','off', 'Position',[scrsz(3)/2-500 scrsz(4)/2-300 1000 400]);
subplot(131);
kx = simplex('OneMExpFun',1,0,[],[],[],1:n,gaussx.^2);
[err, cx, gausstx] = OneMExpFun(kx, 1:n, gaussx.^2);
xlabel('Time');
title('x-coordinate');
ylabel('Width of Gaussian');
ax = axis; ax(4) = ax(4)-ax(3);
text(0.5*n,ax(3) + 0.2*ax(4),['kx = ' num2str(kx)]);
text(0.5*n,ax(3) + 0.1*ax(4),['cx = ' num2str(cx(end))]);
subplot(132);
ky = simplex('OneMExpFun',1,0,[],[],[],1:n,gaussy.^2);
[err, cy, gaussty] = OneMExpFun(ky, 1:n, gaussy.^2);
xlabel('Time');
title('y-coordinate');
ax = axis; ax(4) = ax(4)-ax(3);
text(0.5*n,ax(3) + 0.2*ax(4),['ky = ' num2str(ky)]);
text(0.5*n,ax(3) + 0.1*ax(4),['cy = ' num2str(cy(end))]);


subplot(133);
kz = simplex('OneMExpFun',1,0,[],[],[],1:n,gaussz.^2);
[err, cz, gausstz] = OneMExpFun(kz, 1:n, gaussz.^2);
xlabel('Time');
title('z-coordinate');
ax = axis; ax(4) = ax(4)-ax(3);
text(0.5*n,ax(3) + 0.2*ax(4),['kz = ' num2str(kz)]);
text(0.5*n,ax(3) + 0.1*ax(4),['cz = ' num2str(cz(end))]);

if nargin>1 & ~isempty(flag) & flag==1
   eval(['save ' fname ' n tx tz hx hy hz htx hty htz kx ky kz cx cy cz gaussx gaussy gaussz gausstx gaussty gausstz']);
end
