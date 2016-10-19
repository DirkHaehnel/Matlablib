function [tx, hx, htx, gaussx, gausstx, kx, cx] = FrancescoXLin(x);

% [tx, hx, htx, gaussx, gausstx, kx, cx] = FrancescoXLin(x);
% FrancescoXLin evaluates Borwinian motion data of an optical tweezer
% x - input data
% tx - abscissa of histograms
% hx - histograms
% htx - modeled Gaussians to histograms
% kx - tweezer stiffness
% cx - 4*diffusion/kx

close all
n = 120; % number of binnings
m = 1; % bin width

% determine smalles voltage increment:
tmp = abs(diff(x)); tmp(tmp==0) = [];
dt = min(tmp);

% bulding up the histograms:
dx = x(n*m+1:end)-x(1:end-n*m);
tx = ceil(max(abs(dx))/dt);
tx = dt*(-tx:tx);
hx = [];
for j = 1:n
   dx = x(j*m+1:end) - x(1:end-j*m);
   hx = [hx mhist(dx, tx)];
end

% normalize distributions:
hx = hx./(ones(size(hx,1),1)*sum(hx));
hx = hx(2:end-1,:);
tx = tx(2:end-1);

close all
scrsz = get(0,'ScreenSize'); 
figure('Name','Histogram', 'NumberTitle','off', 'Position',[scrsz(3)/2-500 scrsz(4)/2-100 1000 400]);
plot(tx,hx); axis tight;
xlabel('X');
ylabel('probability');
drawnow

% calculate Gaussian parameters and distributions:
mx = tx*hx;
hh = figure('Name','Gaussian fits', 'NumberTitle','off');
htx = [];
for j = 1:n
	%tmp = simplex('GaussFun',[mx(j) dx(j)],[-inf 0],[],[],[],tx,hx(:,j));
   %gaussx = [gaussx tmp(2)];
   gaussx(j) = sqrt((tx-mx(j)).^2*hx(:,j));
   [tmp, tmp] = GaussFun([mx(j) sqrt(2)*gaussx(j)],tx,hx(:,j));
   htx = [htx tmp];
end
close(hh);

% calculate final parameters:
figure('Name','Parameter fits', 'NumberTitle','off', 'Position',[scrsz(3)/2-500 scrsz(4)/2-300 1000 400]);
kx = simplex('OneMExpFun',1,0,[],[],[],1:n,gaussx.^2);
[err, cx, gausstx] = OneMExpFun(kx, 1:n, gaussx.^2);
xlabel('Time');
title('x-coordinate');
ylabel('Width of Gaussian');
ax = axis; ax(4) = ax(4)-ax(3);
text(0.5*n,ax(3) + 0.2*ax(4),['kx = ' num2str(kx)]);
text(0.5*n,ax(3) + 0.1*ax(4),['cx = ' num2str(cx(end))]);

