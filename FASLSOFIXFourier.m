function res = FASLSOFIXFourier(im,NA,n1,n,n2,d1,d,d2,lambda,mag,pix)

[nx,ny,nz] = size(im);
nn = ceil(0.3*mag/pix);
z = 0;
focpos = 0;
[intx, inty, intz] = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos);
mdf = intx+intx'+inty+inty'+intz+intz';
mdf = [zeros((nx-2*nn-1)/2,size(mdf,2)); mdf; zeros((nx-2*nn-1)/2,size(mdf,2))];
mdf = [zeros(size(mdf,1),(ny-2*nn-1)/2), mdf, zeros(size(mdf,1),(ny-2*nn-1)/2)];
mx = (nx-1)/2;
my = (ny-1)/2;
otm1 = abs(fftshift(fft2(mdf)));
otm2 = mConv2(otm1,otm1);
otm2 = otm2/max(otm2(:));
otm3 = interp2((-2*my:2:2*my)',-2*mx:2:2*mx,otm1,(-my:my)',-mx:mx);
otm3 = otm3/max(otm3(:));

weight = otm3./(0.001+otm2); 
res = im;
for j=1:nz
    res(:,:,j) = abs(fft2(fftshift(ifftshift(ifft2(im(:,:,j))).*weight)));
end

