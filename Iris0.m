function [bck, cx] = Iris(name, nn, bin, conv_flag);

% converting image stacks into movies

if nargin==0
    [filename, pathname] = uigetfile('*.tif', 'Interactive mode data:', 0, 0)
    name = [pathname filename];
end
if nargin<2 | isempty(nn)
    nn = 1;
end
if nargin<3 | isempty(bin)
    bin = 0;
end
imlen = length(imfinfo(name));
mov = moviein(imlen/nn);
close all
tmp = double(imread(name,1));
for j=2:imlen
    tmp = tmp + double(imread(name,j));
end
figure
set(gcf,'units','normalized','position',[0 0 0.9 0.9])
% mim(log(tmp));
% [b,a] = ginput(2);
% a = round(a); b = round(b);
[n,m] = size(tmp);
a = round([0.2 0.8]*n);
b = round([0.2 0.8]*m);
tmp = tmp(a(1):a(2),b(1):b(2));
if bin>0
    bck = mconv2(tmp,disk(bin,1/bin^2))/imlen;
else
    bck = tmp/imlen;
end
cx = [0.8 1.2];
for j=1:nn
    if mod(j-1,nn)==0
        tmp = double(imread(name,j));
    else
        tmp = tmp + double(imread(name,j));
    end
    if mod(j,nn)==0 
        tmp = tmp(a(1):a(2),b(1):b(2));
        if bin>0
            tmp = mconv2(tmp,disk(bin,1/bin^2));
        end
        tmp = tmp./bck;
        mim(tmp);
        caxis(cx)
        mov(:,1) = getframe;        
    end
end
for j=nn+1:imlen
    if mod(j-1,nn)==0
        tmp = double(imread(name,j));
    else
        tmp = tmp + double(imread(name,j));
    end
    if mod(j,nn)==0 
        tmp = tmp(a(1):a(2),b(1):b(2));
        if bin>0
            tmp = mconv2(tmp,disk(bin,1/bin^2));
        end
        tmp = tmp./bck;
        mim(tmp);
        caxis(cx)
        mov(:,j/nn) = getframe;
    end
end
movie2avi(mov,[name(1:end-4) 'bin' num2str(bin) '.avi'],'compression','none');
close all

return

surf(im-abs(ifft2(fftshift(fftshift(fft2(im)).*exp(-0.01*(x.^2+y.^2))))))

