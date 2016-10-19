% HiroshiTest

NA = 1.2;
n0 = 1.52;
n = 1.;
n1 = 1.;
d0 = [];
d = 0; d1 = [];
lamex = 0.62;
mag = 350;
focus = 1;
pixel = 16;

be_res = 10; % minimum resolution of in-plane angle
al_res = 10; % minimum resolution of out-of-plane angle

nn = 20;
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

[intx inty intz, rho, tmp, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 3], 0, NA, n0, n, n1, d0, d, d1, lamex, mag, focus);
for cnt=1:length(theta)
    al = theta(cnt);
    be = -phi(cnt);
    mask(:,:,cnt) = SEPImage(al,be,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
    cnt
end

signal = 5e3;
bckgrnd = 0;
for cnt=1:length(theta)
    for k=1:10
        int = mask(:,:,cnt);
        int = int/sum(int(:));
        z = poissrnd(bckgrnd+signal*int);
        imwrite(uint16(z),'HiroshiTestFile.tiff','compression','none','writemode','append');
    end
end

return

res = Hiroshi('HiroshiTestFile.tiff',0,[1 1; 2*nn+1 2*nn+1]);
