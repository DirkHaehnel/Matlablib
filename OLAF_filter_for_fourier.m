% makes a filter for Fourier reweighting, cutting off all "impossible"
% frequencies

%requires a (square) image file (im) with an even number of pixels along both
%directions


lamex = 625;  % excitation wavelength in nm (for widefield put same as emission)
lamem = 625;  % emission wavelength in nm 
pixel = 50;  %pixel size in nm
NA = 1.4;   %NA of the objective

size = (size(im,1));
lambda = (lamex+lamem)/2;
k_max = 4*pi*NA*pixel/lambda;
fourierpixel = 2*pi/size;
k_max_pixel = k_max/fourierpixel;

c = k_max_pixel+20;
a = Disk(c);
b = zeros(size);
e = size/2;
b(e-c:e+c , e-c:e+c) = a;

[x,y] = meshgrid(-20:20);
f = exp(-(x.^2+y.^2)/10^2/2);
d = conv2(b,f,'same');
filter = d/max(max(d));


