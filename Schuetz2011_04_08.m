%TIRF imaging PSF

close all
clear all

pix = 16;
nn = 15;
zv = (0:10)*50/1e3;
NA = 1.45;
n1 = 1.52;
n = 1.33;
n2 = 1.33;
d1 = [];
d = 0;
d2 = [];
lambda = 0.55;
mag = 150;
atf = [];
ring = [];

for j=1:length(zv)
    z = zv(j);
    focpos = z;
    [intx, inty, intz] = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring);
    int(:,:,j) = intx+intx'+inty+inty'+intz+intz';
end
CombineImages(int,2,5,[],{'0 nm','50 nm','100 nm','150 nm','200'},{'+ 0 nm','+ 250 nm'})

