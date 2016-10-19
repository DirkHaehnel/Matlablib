function [im1, im2, amp1, amp2, xc, yc, ww, im] = SarahFRET(name, ww)

tic

close all 

if nargin==0
    [filename, pathname] = uigetfile('*.tif', 'Interactive mode data:', 0, 0);
    name = [pathname filename];
end

fac = 30;
tsh = 1;

% read data
close all
tmp = imfinfo(name);
imlen = length(tmp);
im = zeros(tmp(1).Height,tmp(1).Width,imlen);
for j=1:imlen
    tmp = imread(name,j);
    im(:,:,j) = double(tmp);
end

toc

% split image
[n,m,imlen] = size(im);
if mod(m,2)==1
    im1 = im(:,1:(end-1)/2,:);
    im2 = im(:,(end-1)/2+1:end-1,:);
    m = (m-1)/2;
else
    im1 = im(:,1:end/2,:);
    im2 = im(:,end/2+1:end,:);
    m = m/2;
end

if 1
    tic    
    
    % background reduction
    [n,m,imlen] = size(im1);
    [x,y] = meshgrid(-(m-1)/2:(m-1)/2,-(n-1)/2:(n-1)/2);
    for j=1:imlen
        bck1 = abs(ifft2(ifftshift(fftshift(fft2(im1(:,:,j))).*exp(-(x.^2+y.^2)/fac))));
        im1(:,:,j) = im1(:,:,j) - bck1;
        bck2 = abs(ifft2(ifftshift(fftshift(fft2(im2(:,:,j))).*exp(-(x.^2+y.^2)/fac))));
        im2(:,:,j) = im2(:,:,j) - bck2;
    end
    
    toc
    
    tic
    
    % cut boundaries
    boundary_width = 5;
    a = round([boundary_width n-boundary_width]);
    b = round([boundary_width m-boundary_width]);
    im1 = im1(a(1):a(2),b(1):b(2),:);
    im2 = im2(a(1):a(2),b(1):b(2),:);
    bck1 = bck1(a(1):a(2),b(1):b(2));
    bck2 = bck2(a(1):a(2),b(1):b(2));
    [n,m,imlen] = size(im1);
        
    toc
    
end

[x,y] = meshgrid(-(m-1)/2:(m-1)/2,-(n-1)/2:(n-1)/2);

tic

% image registration
mm = 8;
imm1 = mean(im1,3); imm2 = mean(im2,3);
for j=1:2*mm+1
    for k=1:2*mm+1
        tst(j,k) = sum(sum(imm1(mm+1:end-mm,mm+1:end-mm).*imm2(j:n-2*mm-1+j, k:m-2*mm-1+k)));
    end
end
[xx,yy] = meshgrid(-mm:mm,-mm:mm);
ind = tst==max(tst(:));
for j=1:imlen
    im2(:,:,j) = interp2(x,y,im2(:,:,j),x+xx(ind),y+yy(ind),'cubic',0);
end

% cut boundaries
boundary_width = 5;
a = round([boundary_width n-boundary_width]);
b = round([boundary_width m-boundary_width]);
im1 = im1(a(1):a(2),b(1):b(2),:);
im2 = im2(a(1):a(2),b(1):b(2),:);
[n,m,imlen] = size(im1);
[x,y] = meshgrid(-(m-1)/2:(m-1)/2,-(n-1)/2:(n-1)/2);

toc

% determine size of single molecule image
if nargin<2 || isempty(ww)
    im = mean(im1+im2,3);
    mim(im)
    [b,a] = ginput(2);
    a = round(a);b=round(b);
    [mx,my,wx,wy] = GaussEllipse([],[],abs(im(a(1):a(2),b(1):b(2))));
    ww = sqrt(wx.*wy)/2;
end

% find molecule positions
w3 = round(3*ww);
[xx,yy] = meshgrid(-w3:w3,-w3:w3);
mask = exp(-(xx.^2+yy.^2)/ww^2/2);
[tmp, tmp, tmp, tmp, xc, yc] = FindPattern(mean(im1+im2-min(im1(:)+im2(:)),3), mask, [], [], 1, tsh); 
subplot(121); mim(mean(im1+im2,3)); hold on; plot(xc,yc,'oc'); hold off; axis on
subplot(122); imagesc(cat(3,(mean(im1,3)-min(im1(:)))/(max(im1(:))-min(im1(:))),(mean(im2,3)-min(im2(:)))/(max(im2(:))-min(im2(:))),zeros(size(im1,1),size(im1,2)))); axis image

% determine brightness values
amp1 = zeros(imlen,length(xc));
amp2 = amp1;
for j=1:imlen
    for k=1:length(xc)
        amp1(j,k) = sum(sum(im1(yc(k)-w3:yc(k)+w3,xc(k)-w3:xc(k)+w3,j)));
        amp2(j,k) = sum(sum(im2(yc(k)-w3:yc(k)+w3,xc(k)-w3:xc(k)+w3,j)));
    end
end

% sort data
[ind,ind] = sort(sum(amp1+amp2));
amp1 = amp1(:,ind);
amp2 = amp2(:,ind);
xc = xc(ind);
yc = yc(ind);
ind = sum(amp1+amp2)>0;
amp1 = amp1(:,ind);
amp2 = amp2(:,ind);
xc = xc(ind);
yc = yc(ind);

% show final result
figure
mim([amp1; inf*ones(1,size(amp1,2)); amp2]);
colormap jet

return

% Anzeige der Intensitäten
t = 1:size(amp1,1); 
for j=1:size(amp1,2) plot(t,amp1(:,j),t,amp2(:,j)); title(int2str(j)); pause; end

% Anzeigen Position des z.B. 54ten Moleküls:
k=54; subplot(121); mim(mean(im1,3)); hold on; plot(xc(k),yc(k),'oc'); hold off; subplot(122); mim(mean(im2,3)); hold on; plot(xc(k),yc(k),'oc'); hold off

% Alle Moleküle anzeigen:
subplot(121); mim(mean(im1,3)); hold on; plot(xc,yc,'oc'); hold off; subplot(122); mim(mean(im2,3)); hold on; plot(xc,yc,'oc'); hold off

% Colocalization plot
imagesc(cat(3,(mean(im1,3)-min(im1(:)))/(max(im1(:))-min(im1(:))),(mean(im2,3)-min(im2(:)))/(max(im2(:))-min(im2(:))),zeros(size(im1,1),size(im1,2)))); axis image
