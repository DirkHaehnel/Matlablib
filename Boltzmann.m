function [tx, ty, tw, tz, hx, hy, hw, hz, cx, cy] = Boltzmann(fname, flag);

% [tx, ty, tw, tz, hx, hy, hw, hz, cx, cy] = Boltzmann(fname)
% Boltzmann calculates the position histograms of the Borwinian motion data of an optical tweezer
% fname - name of input data file
% t* - axes of histogramming
% h* - histograms
% c* - quadratic potential coefficients

close all
intake = 5e4; % how many data to process at once

if nargin<2
   flag = 0;
else
   flag = 1;
end

fin = fopen(fname, 'r');
[data, count] = fscanf(fin, '%f ', intake);

an = max(union(0,union(findstr(fname,'\'),findstr(fname,'/'))))+1;
en = min(union(findstr(fname,'.'),length(fname)+1))-1;
fname = fname(an:en) ;

% determine smalles voltage increment:
tmp = abs(diff(data));
tmp(tmp==0) = [];
dt = min(tmp);

% bulding up the histograms:
x = data(1:5:end) - data(2:5:end) - data(3:5:end) + data(4:5:end); 
y = data(1:5:end) + data(2:5:end) - data(3:5:end) - data(4:5:end);
w = abs(data(1:5:end) + data(2:5:end) + data(3:5:end) + data(4:5:end));
tx = ceil(max(abs(x))/dt);
ty = ceil(max(abs(y))/dt);
z = data(5:5:end);
if flag
	x = x./w;
	y = y./w;
	tx = dt*(-tx:tx)/mean(w);
	ty = dt*(-ty:ty)/mean(w);
else
   tx = dt*(-tx:tx);
	ty = dt*(-ty:ty);
end
tz = dt*(floor(min(z)/dt):ceil(max(z)/dt));
tw = dt*(floor(min(w)/dt):ceil(max(w)/dt));

hx = mhist(x, tx);
hy = mhist(y, ty);
hw = mhist(w, tw);
hz = mhist(z, tz);
scrsz = get(0,'ScreenSize'); 
figure('Name',['Histograms ' fname], 'NumberTitle','off', 'Position',[scrsz(3)/2-500 scrsz(4)/2-100 1000 400]);
subplot(121);
plot(tx,hx./(ones(size(hx,1),1)*sum(hx)),ty,hy./(ones(size(hy,1),1)*sum(hy))); axis tight;
xlabel('X/Y');
ylabel('probability');
subplot(122);
plot(tz,hz./(ones(size(hz,1),1)*sum(hz))); axis tight;
xlabel('Z');
drawnow
while count>0
   [data, count] = fscanf(fin, '%f ', intake);
   if count>0
      x = data(1:5:end) - data(2:5:end) - data(3:5:end) + data(4:5:end);
      y = data(1:5:end) + data(2:5:end) - data(3:5:end) - data(4:5:end);
		w = abs(data(1:5:end) + data(2:5:end) + data(3:5:end) + data(4:5:end));
      z = data(5:5:end);
      if flag
         x = x./w; 
         y = y./w;
      end
      hx = hx + mhist(x, tx);
      hy = hy + mhist(y, ty);
      hw = hw + mhist(w, tw);
      hz = hz + mhist(z, tz);
      subplot(121);
      plot(tx,hx./(ones(size(hx,1),1)*sum(hx)),ty,hy./(ones(size(hy,1),1)*sum(hy))); axis tight;
      xlabel('X/Y');
      ylabel('probability');
      subplot(122);
      plot(tz,hz./(ones(size(hz,1),1)*sum(hz))); axis tight;
      xlabel('Z');
      drawnow
   end
end
fclose(fin);

% normalize distributions:
hx = hx./(ones(size(hx,1),1)*sum(hx));
hy = hy./(ones(size(hy,1),1)*sum(hy));
hw = hw./(ones(size(hw,1),1)*sum(hw));
hz = hz./(ones(size(hz,1),1)*sum(hz));
hx = hx(2:end-1);
hy = hy(2:end-1);
hw = hw(2:end-1);
hz = hz(2:end-1);
tx = tx(2:end-1);
ty = ty(2:end-1);
tw = tw(2:end-1);
tz = tz(2:end-1);
tx(hx==0) = [];
ty(hy==0) = [];
tz(hz==0) = [];
tw(hw==0) = [];
hx(hx==0) = [];
hy(hy==0) = [];
hz(hz==0) = [];
hw(hw==0) = [];

% plot final results
cx = polyfit(tx,-log10(hx'),2);
cy = polyfit(ty,-log10(hy'),2);
subplot(121)
plot(tx,-log10(hx),tx,polyval(cx,tx)); axis tight;
xlabel('X');
ylabel('potential');
subplot(122)
plot(ty,-log10(hy),ty,polyval(cy,ty)); axis tight;
xlabel('Y');
set(gcf,'Name',['Potentials ' fname]);
drawnow
