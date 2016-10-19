function [tx, hx, htx, delay, gaussx, gausstx, kx, cx] = FrancescoX(x,s);

% [tx, hx, htx, delay, gaussx, gausstx, kx, cx] = FrancescoX(x,s);
% FrancescoX evaluates Borwinian motion data of an optical tweezer
% x - input data
% s - sum signal
% tx - abscissa of histograms
% hx - histograms
% htx - modeled Gaussians to histograms
% delay - delay time values
% gaussx = Gaussian width values
% gausstx - fitted Gaussian width values
% kx - tweezer stiffness
% cx - 4*diffusion/kx

close all
m = 8; % number of subdivions
n = 4; % number of cascades

% determine smalles voltage increment:
tmp = abs(diff(x));
tmp(tmp==0) = [];

if nargin>1
    dt = min(tmp)/mean(s);
    x = x./s;
else
    dt = min(tmp);
end

% bulding up the histograms:
dx = x(m^(n-1)+1:end)-x(1:end-m^(n-1));
tx = ceil(max(abs(dx))/dt);
tx = dt*(-tx:tx);
hx = []; delay = [];
for j = 0:n-1
    for k = 1:m-1
        delay = [delay k*m^j];
        dx = x(k*m^j+1:end) - x(1:end-k*m^j);
        hx = [hx mhist(dx, tx)];
    end
end
close all
scrsz = get(0,'ScreenSize'); 
figure('Name','Histograms', 'NumberTitle','off', 'Position',[scrsz(3)/2-500 scrsz(4)/2-100 1000 400]);
plot(tx,hx./(ones(size(hx,1),1)*sum(hx))); axis tight;
xlabel('X');
ylabel('probability');
drawnow

% normalize distributions:
hx = hx./(ones(size(hx,1),1)*sum(hx));
hx = hx(2:end-1,:);
tx = tx(2:end-1);

% plot final results
plot(tx,hx);
xlabel('X');
ylabel('probability');
drawnow

mx = tx*hx;
hh = figure('Name','Gaussian fits', 'NumberTitle','off');
htx = [];
for j = 1:length(delay)
   gaussx(j) = sqrt((tx-mx(j)).^2*hx(:,j));
   [tmp, tmp] = GaussFun([mx(j) sqrt(2)*gaussx(j)],tx,hx(:,j));
   htx = [htx tmp];
end
close(hh);

% calculate final parameters:
figure('Name','Parameter fits', 'NumberTitle','off', 'Position',[scrsz(3)/2-500 scrsz(4)/2-300 1000 400]);
kx = simplex('OneMExpFun',1,0,[],[],[], delay, gaussx.^2);
[err, cx, gausstx] = OneMExpFun(kx, delay, gaussx.^2);
xlabel('Time');
title('x-coordinate');
ylabel('Width of Gaussian');
ax = axis; ax(4) = ax(4)-ax(3);
text(0.5*delay(end),ax(3) + 0.2*ax(4),['kx = ' num2str(kx)]);
text(0.5*delay(end),ax(3) + 0.1*ax(4),['cx = ' num2str(cx(end))]);


