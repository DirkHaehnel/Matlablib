function [res, tmp1, tmp2] = FourierUpsampling(im,nraster,flag)

if nargin<3 || isempty(flag)
    m = size(im,1);
    n = size(im,2);
    tmp1 = [im(floor(m/2):-1:1,floor(n/2):-1:1) im(floor(m/2):-1:1,:) im(floor(m/2):-1:1,end:-1:end-ceil(n/2)+1); ...
        im(:,floor(n/2):-1:1) im im(:,end:-1:end-ceil(n/2)+1); ...
        im(end:-1:end-ceil(m/2)+1,floor(n/2):-1:1) im(end:-1:end-ceil(m/2)+1,:) im(end:-1:end-ceil(m/2)+1,end:-1:end-ceil(n/2)+1)];
else
    tmp1 = im;
end

fac = nraster;
a = size(tmp1,1);
b = size(tmp1,2);
nraster = [ceil(ceil((a-1)/2)*(nraster-1)) ceil(ceil((b-1)/2)*(nraster-1))];

tmp2 = real(ifft2(ifftshift(padarray(fftshift(fft2(tmp1)),nraster))));
tmp2 = circshift(tmp2,ceil([fac,fac]/2-1));

if nargin<3 || isempty(flag)
    res = tmp2(floor(size(tmp2,1)/4)+1:end-ceil(size(tmp2,1)/4),floor(size(tmp2,2)/4)+1:end-ceil(size(tmp2,2)/4));
else
    res = tmp2;
end

return

[x,y]=meshgrid(1:7,1:7);
im=cos(0.5*(y-x));
mim(im)
[res1,tmp11,tmp12]=FourierUpsampling(im,2,1);
[res2,tmp21,tmp22]=FourierUpsampling(im,2);
