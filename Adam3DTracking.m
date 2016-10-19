% program Adam3DTracking

close all

if ~exist('x')==1
    pth = 'c:\Joerg\Doc\Microscopy\Tracking3D\aug3 qd1\';
    filenames = dir([pth '*.tif']);
    for j=1:length(filenames)
        x(:,:,j) = imread([pth filenames(j).name]);
    end
end

NA = 1.0;
n0 = 1.51;
n = 1.4;
n1 = 1.4;
d0 = [];
d = 0;
d1 = [];
lamem = 0.57;
mag = 600;
pixel = 16;
pic = 0;
be_res = [];
al_res = [];
nn = 499.5;
focpos = 0;
zv = 0:0.1:20;

rad = (nn+1)*sqrt(2)*pixel/mag;
rho = (0.5:1000)/1000*rad;
int = zeros((2*nn+1)/10,(2*nn+1)/10,length(focusv));
for j=1:length(focusv)
    z = zv(j);
    [intx inty intz, rho, tmp, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 1.5*max(nn)*pixel/mag], z, NA, n0, n, n1, d0, d, d1, lamem, mag, focpos);
    mask = SEPImage(pi/2,0,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
    mask = mask + SEPImage(pi/2,pi/2,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
    mask = mask + SEPImage(0,0,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
    for kx=1:(2*nn+1)
        for ky = 1:(2*nn+1)
            int(ceil(kx/10),ceil(ky/10),j) = int(ceil(kx/10),ceil(ky/10),j) + mask(kx,ky);
        end
    end
    int(:,:,j) = int(:,:,j)/sum(sum(int(:,:,j)));
    surf(-nn/10:nn/10,-nn/10:nn/10,int(:,:,j)); drawnow
end


return

for k=1:size(x,3) 
    [err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm(:,:,k)] = FindPattern(double(x(:,:,k)),int,[],[],[], 1); 
end
