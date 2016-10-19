NA = 1.45;
n0 = 1.52;
n = 1.333;
n1 = 1.333;
d0 = [];
d = 0; d1 = [];
lamex = 0.575;
mag = 160;
focus = 0.5;

nn = 20; % 2*nn+1 is the edge size of the generated frames
n2 = 2*nn;
n3 = (2*nn+1)^2;
[x,y] = meshgrid(-nn:nn,-nn:nn);
p = angle(x+i*y);
r = sqrt(x.^2+y.^2);

[intx inty intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 3], 0, NA, n0, n, n1, d0, d, d1, lamex, mag, focus);
rho = rho/16;

al_step = 10; % step size for polar angle in degrees
be_step = 10; % step size for azimuthal angle in degrees
signal = 5e3; % average total number of detected photons per pattern
bckgrnd = 0; % how many background photons
max_cnt = 1000; % max number of frames per sequence

for j=1:round(90/al_step)
    al = (j-1)*al_step/180*pi;
    if (j-1)*al_step==90 kmax=1; else kmax = round(90/be_step); end
    for k=1:kmax
        be = (k-1)*be_step/180*pi;
        outfile = ['TestIm_alpha_' mint2str((j-1)*al_step,2) '_beta_' mint2str((k-1)*be_step,3) '.tif'];
        for jj=1:max_cnt
            x0 = rand-0.5; y0 = rand-0.5;
            r = sqrt((x-x0).^2+(y-y0).^2);
            p = angle((x-x0)+i*(y-y0));
            int = real((cos(al)*(interp1(rho,fxx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,fxx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,fxz,r,'cubic')).*...
                conj(cos(al)*(interp1(rho,byx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,byx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,byz,r,'cubic')) + ...
                (cos(al)*sin(2*(p-be)).*interp1(rho,fxx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,fxz,r,'cubic')).*...
                conj(cos(al)*sin(2*(p-be)).*interp1(rho,byx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,byz,r,'cubic')));
            int = int/sum(int(:));
            z = poissrnd(bckgrnd+signal*int);
            z(isnan(z)) = 0;
            imwrite(uint16(z),outfile,'tif','compression','none','writemode','append');
        end
    end
end

