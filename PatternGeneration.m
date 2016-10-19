function model = PatternGeneration(z, NA, n0, n, n1, d0, d, d1, lamem, mag, focus, atf, ring, pixel, nn, be_res, al_res, pic)
                 
if nargin<11
    atf = [];
end
if nargin<12
    ring = [];
end
if nargin<18 || isempty(pic)
    pic = 0;
end
if nargin<16 || isempty(be_res)
    be_res = 10; % minimum resolution of in-plane angle
end
if nargin<17 || isempty(al_res)
    al_res = 10; % minimum resolution of out-of-plane angle
end
if nargin<15 || isempty(nn)
    nn = [10 10];
end
bck = Disk(nn);
cnt = 1;
for k=90:-al_res:0
    al = k/180*pi;
    if k==90
        jj = round(180/be_res);
        dbe = pi/jj;
    elseif k==0
        jj = 1;
        dbe = 0;
    else
        jj = round(sin(al)*360/be_res);
        dbe = 2*pi/jj;
    end
    for j=1:jj
        theta(cnt) = al;
        phi(cnt) = dbe*(j-1);
        cnt = cnt+1;
    end
end

[intx inty intz, rho, tmp, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 1.5*max(nn)*pixel/mag], z, NA, n0, n, n1, d0, d, d1, lamem, mag, focus, atf, ring);
for cnt=1:length(theta)
    al = theta(cnt);
    be = -phi(cnt);
    mask(:,:,cnt) = SEPImage(al,be,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
end
if pic==1
    col = ceil(sqrt(size(mask,3)));
    wdth = size(mask,2);
    hght = size(mask,1);
    im = zeros(ceil(size(mask,3)/col)*hght,col*wdth);
    for j=1:size(mask,3)
        im(fix((j-1)/col)*hght+1:(fix((j-1)/col)+1)*hght,mod(j-1,col)*wdth+1:(mod(j-1,col)+1)*wdth) = mask(:,:,j);
    end
    mim(im)
end

for j=1:size(mask,3) 
    mask(:,:,j) = mask(:,:,j).*bck; 
    mask(:,:,j) = mask(:,:,j)/sum(sum(mask(:,:,j))); 
end

model.rho = rho;
model.theta = theta;
model.phi = phi;
model.mask = mask;
model.fxx0 = fxx0;
model.fxx2 = fxx2;
model.fxz = fxz;
model.byx0 = byx0;
model.byx2 = byx2;
model.byz = byz;
