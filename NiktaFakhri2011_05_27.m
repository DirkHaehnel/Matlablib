fname = 'c:\Joerg\Doc\Microscopy\NiktaFakhri\Microrheology\NTinactin.tif';
imnum = length(imfinfo(fname));

im = double(imread(fname,1));
for j=2:imnum
    im = im + double(imread(fname,j));
end

[x,y]=meshgrid(-size(im,2)/2:(size(im,2)/2-1),-size(im,1)/2:(size(im,1)/2-1));
imf = fftshift(fft2(im-mean(im(:))));

indx = ~x(1,:)==0;
indy = ~y(:,1)==0;
surf(log(imf(indy,indx)))

r = sqrt(x.^2+y.^2);

if 0
    z = besselj(1,r)./r;
    z(isnan(z)) = 1/2;
    z = abs(fftshift(fft2(abs(z))));
    p = Simplex('AffineFit2D', [0 0 0.01], [0 0 0], [0 0 inf], [], [], x(indy,indx), y(indy,indx), abs(imf(indy,indx)), x(indy,indx), y(indy,indx), z(indy,indx), 1)
end

if 1
    z = exp(-r.^2/2/30^2);
    bck = real(ifft2(ifftshift(imf.*(1-z))));
    imm = real(ifft2(ifftshift(imf.*z)));
    CombineImages(cat(3,im,imm-bck),1,2)
end
    



