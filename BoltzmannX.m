function [tx, hx, cx] = BoltzmannX(x,s);

% [tx, hx, cx] = BoltzmannX(x,s)
% Boltzmann calculates the position histograms of the Borwinian motion data of an optical tweezer
% x - difference signal
% s - sum signal
% tx - axis of histogramming
% hx - histogram
% cx - quadratic potential coefficients

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
tx = ceil(max(abs(x))/dt);
tx = dt*(-tx:tx);

hx = mhist(x, tx)';
scrsz = get(0,'ScreenSize'); 
figure('Name','Histogram', 'NumberTitle','off', 'Position',[scrsz(3)/2-500 scrsz(4)/2-100 1000 400]);
plot(tx,hx./(ones(size(hx,1),1)*sum(hx))); axis tight;
xlabel('X');
ylabel('probability');
drawnow

% normalize distributions:
hx = hx/sum(hx);
hx = hx(2:end-1);
tx = tx(2:end-1);
tx(hx==0) = [];
hx(hx==0) = [];

% plot final results
cx = polyfit(tx,-log10(hx),2);
plot(tx,-log10(hx),tx,polyval(cx,tx)); axis tight;
xlabel('X');
ylabel('potential');
set(gcf,'Name','Potential');
drawnow
