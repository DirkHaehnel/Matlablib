bckgrnd = 0; % number of background photons per pixel
signal = 1e4; % total average number of fluorescence photons 
pixel = 43; % pixel size in nm
max_cnt = 1e4; % maximum number of generated frames
al = pi/4; % incliniation angle of molecule dipole

NA = 1.45; % NA of objective
n0 = 1.52;
n = 1.333;
n1 = 1.333;
d0 = [];
d = 0; d1 = [];
lamex = 0.575; % emission wavelength
mag = 16e3/pixel; % 16 mum / pixel size in nm
focus = 0; % defocusing in mum

nn = 15; 
n2 = 2*nn+1; % frame size
[x,y] = meshgrid(-nn:nn,-nn:nn);
p = angle(x+i*y);
r = sqrt(x.^2+y.^2);

[intx inty intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 3], 0, NA, n0, n, n1, d0, d, d1, lamex, mag, focus);
rho = rho/16;
be = 0;
int = real((cos(al)*(interp1(rho,fxx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,fxx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,fxz,r,'cubic')).*...
    conj(cos(al)*(interp1(rho,byx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,byx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,byz,r,'cubic')) + ...
    (cos(al)*sin(2*(p-be)).*interp1(rho,fxx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,fxz,r,'cubic')).*...
    conj(cos(al)*sin(2*(p-be)).*interp1(rho,byx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,byz,r,'cubic')));    
int = int/sum(int(:));
z = zeros(n2,n2,max_cnt);
for j=1:max_cnt
    z(:,:,j) = poissrnd(bckgrnd+signal*int);
    if mod(j,100)==0 j, end
end



