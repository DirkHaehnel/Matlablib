% program for calculating defoc excitation scan images of a dipole emitter
% within a diamond film
%
% all length values are in units of micrometers!!

lamex = 0.76; % excitation wavelength
resolution = 50;
rhofield = [-lamex/resolution(1)/2 3.];
molpos = 4; % molecule's distance in diamond layer from surface
pixel = 0.01; % CCD pixel size
NA = 1.3; % N.A. of objective
fd = 3e3; % foccal distance of objective (Olympus)
n0 = 1.51; % immersion oil ref. index
n = 2.417; % diamond ref. index
n1 = n;
d0 = [];
d = 0;
d1 = [];
over = inf; % diffraction limited focusing
focposv = molpos*n0/n - (0:47)*0.05; % vector of defoc values (should be 48 values!!, currently in steps of 50 nm)
atf = [];
% out of plane rotation of dipole axis
theta = (90:-5:40)/180*pi;
% in-plane rotation of dipole axis
phi = zeros(size(theta)); 
psi = (180-24)/180*pi;

nn = 100;
mask = zeros(2*nn+1,2*nn+1,length(focposv),length(theta));
for jz=1:length(focposv)
    [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, molpos, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposv(jz), atf, resolution);
    for cnt=1:length(theta)
        al = theta(cnt);
        be = phi(cnt);
        mask(:,:,jz,cnt) = DefExcImage(al,be,nn,pixel,rho,fxc,fxs,fyc,fys,fzc,fzs,psi);
    end
end
for j=1:8 sx{j} = [int2str((max(focposv)-focposv(j))*1e3) ' nm']; end
for j=1:6 sy{j} = ['+' int2str((focposv(1)-focposv(9))*1e3) ' nm']; end
sy{1} = '';

% shows the defoc images for different defoc values for an inclinatin angle
% of theta(11); change '11' to other values to see result for other inclination angles
CombineImages(mask(:,:,:,11),6,8,'scale',sx,sy);

