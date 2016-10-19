% simulation of a brightness imaging setup

average_photons = 1e3;
nn = 20;
z = -2:0.05:2;
mag = 100; % 
pixel = 6.45;
NA = 1.2;
lamem = 0.67;
n0 = 1.33;
n = n0;
n1 = n0;
d0 = [];
d = 0;
d1 = [];

int = zeros(2*nn+1,2*nn+1,length(z));
for j=1:length(z) 
    int(:,:,j) = DefocImage(0, 0, 0, 0, 1/2, z(j), mag, NA, lamem, n0, n, n1, d0, d, d1); 
    int(:,:,j) = int(:,:,j)/sum(sum(int(:,:,j)));
end

k = 18;
[x,y,wx,wy,amp] = Gauss2D((-nn:nn)*pixel/mag, (-nn:nn)*pixel/mag, poissrnd(int(:,:,k)*average_photons));

return

% movie
for j=1:1000 mim(poissrnd(int(:,:,11)*average_photons)); end