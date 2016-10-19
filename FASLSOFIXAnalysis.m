function [sof, im0] = FASLSOFIXAnalysis(im,ncum,win)

if nargin<2 || isempty(ncum)
    ncum = 4;
end

if nargin<3 || isempty(win)
    win = 1e2;
end

if ncum>6
    error('SOFIAnalysis:argChk', 'Higher cumulants than 4th order are stupid')
end

if ischar(im)
    x = tiffread(im);
    [m,n] = size(x(1).data);
    im = zeros(m,n,length(x));
    for j=1:length(x)
        im(:,:,j) = double(x(j).data); 
        disp(j) 
    end
    clear x
end

% im = im - repmat(mean(im,3),[1 1 size(im,3)]);
% for j=1:ncum-1
%     soffull{j} = zeros(2^j*(size(im,1)-1)+1,2^j*(size(im,2)-1)+1,floor(size(im,3)/win));
% end
[x,y] = meshgrid(0:size(im,2)-1,0:size(im,1)-1);
[xi,yi] = meshgrid(0:1/ncum:(size(im,2)-1),0:1/ncum:(size(im,1)-1));
sof = zeros(ncum*(size(im,1)-1)+1,ncum*(size(im,2)-1)+1,ncum-1);
im0 = zeros(ncum*(size(im,1)-1)+1,ncum*(size(im,2)-1)+1);
tmp = zeros(ncum*(size(im,1)-1)+1,ncum*(size(im,2)-1)+1,win);
for k=1:floor(size(im,3)/win)
    for jj=1:win
        tmp(:,:,jj) = interp2(x,y,im(:,:,(k-1)*win+jj),xi,yi,'spline');
    end
    im0 = im0 + sum(tmp,3); 
    tmp = tmp - repmat(mean(tmp,3),[1 1 size(tmp,3)]);
    for j=1:ncum-1
        sof(:,:,j) = sof(:,:,j) + FASLcumulant0(tmp,j+1,1);
        disp([k j]);
    end
end

mim(cat(3,im0,abs(sof(:,:,1))))
colormap jet
% mim(cat(4,cat(3,sum(im0,3),abs(sof(:,:,1))),abs(sof(:,:,2:3))))

return

x = tiffread('c:\Joerg\Doc\Microscopy\SOFI\PicoQuant\Alexa532_Mikrotubuli_Kammer2_3\Alexa532_Mikrotubuli_Kammer2_3_raw_data.tif');
y = zeros(205,131,1e4);
for j=1:1e4 y(:,:,j)=double(x(j).data); j, end
clear x

mim([(im1-min(min(im1)))/(max(max(im1))-min(min(im1))) sum(sof1,3)/max(max(sum(sof1,3)))]); colormap jet

return

x=-100:100;
z=1-0.1*abs(x); z(z<0)=0;
zz=conv(z,z);
zz=zz(length(z)/2:end-length(z)/2);
y2=real(fftshift(ifft(ifftshift(zz))));
y1=real(fftshift(ifft(ifftshift(z))));
y1=y1/max(y1);
y2=y2/max(y2);
b=ones(size(x)); b(z==0)=0;
bb=conv(b,b);
bb=bb(length(z)/2:end-length(z)/2);
b2=real(fftshift(ifft(ifftshift(bb))));
b1=real(fftshift(ifft(ifftshift(b))));
b1=b1/max(b1);
b2=b2/max(b2);
plot(x,y2,x,b2,'o')
[sqrt(sum(x.^2.*y2)/sum(y2)) sqrt(sum(x.^2.*b2)/sum(b2))]