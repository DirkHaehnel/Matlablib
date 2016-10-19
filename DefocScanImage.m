lamex = 0.64;
resolution = 50;
rhofield = [-lamex/resolution(1)/2 3.];
molpos = 0;
pixel = 0.01; %/sqrt(2);
NA = 1.49;
fd = 3e3;
n0 = 1.51;
n = 1.0;
n1 = 1.0;
d0 = [];
d = 0;
d1 = [];
over = inf;
focposv = (0:47)*0.02;
atf = [];
ring = [];
psi = 0;

theta = (90:-5:0)/180*pi;
% in-plane rotation of dipole axis
phi = 0*pi/2; %(0:5:355)/180*pi;
theta = ones(length(phi),1)*theta;
theta = theta(:)';
phi = repmat(phi,[1 length(theta)/length(phi)]);
nn = 100; %floor(rhofield(2)/pixel/sqrt(2));
mask = zeros(2*nn+1,2*nn+1,length(focposv),length(theta));
for jz=1:length(focposv)
    [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, molpos, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposv(jz), atf, resolution);
    
    for cnt=1:length(theta)
        al = theta(cnt);
        be = phi(cnt);
        mask(:,:,jz,cnt) = DefExcImage(al,be,nn,pixel,rho,fxc,fxs,fyc,fys,fzc,fzs,psi);
    end
end

for j=1:8 sx{j} = [int2str(focposv(j)*1e3) ' nm']; end
if focposv(9)>0
    for j=1:6 sy{j} = ['+' int2str(focposv(9)*1e3) ' nm']; end
else
    for j=1:6 sy{j} = ['-' int2str(abs(focposv(9)*1e3)) ' nm']; end
end
sy{1} = '';
CombineImages(mask(:,:,:,11),6,8,'scale',sx,sy);

