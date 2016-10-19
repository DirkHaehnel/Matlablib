z = 0;
NAv = [1.3 1.4 1.45];
nsio = 1.46;
n = 1.33;
n2 = 1.33;
dag = 30e-3;
dsio = (5:5:50)*1e-3;
d = 0;
d2 = [];
lambda = 0.67;
mag = 300;
focpos = 0.25;

load metals
nag = silver(wavelength==lambda*1e3);
n1 = [1.52 nag nsio];

nn = 20;
pix = 16;

cnt = 1;
intx = zeros(2*nn+1,2*nn+1,numel(NAv)*numel(dsio));
intz = intx;
for j=1:numel(NAv)
    NA = NAv(j);
    for k=1:numel(dsio)
        d1 = [dag dsio(k)];
        [intx(:,:,cnt), ~, intz(:,:,cnt)] = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos); 
        cnt = cnt+1;
    end
end

for j=1:numel(NAv) 
    sx{j}=[mnum2str(NAv(j),1,2)]; 
end; 
for j=1:numel(dsio) 
    sy{j}=[num2str(1e3*dsio(j)) ' nm']; 
end;
CombineImages(intx,numel(NAv),numel(dsio),[],sy,sx)
title({'horizontal orientation',''})
print -dpng -r300 AlexeySilverDefocX

CombineImages(intz,numel(NAv),numel(dsio),[],sy,sx)
title({'vertical orientation',''})
print -dpng -r300 AlexeySilverDefocZ


focposv = [0 0.25 0.5 0.75 1];
theta = (0:5:90)/180*pi;
NA = 1.45;
d1 = [dag dsio(1)];
cnt = 1;
int = zeros(2*nn+1,2*nn+1,numel(focposv)*numel(theta));
for j=1:numel(focposv)
    focpos = focposv(j);
    for k=1:numel(theta)
        int(:,:,cnt) = DefocImage(nn, pix, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, [], [], [theta(k) 0]);
        cnt = cnt+1;
    end
end
for j=1:numel(theta) s{j} = [num2str(theta(j)/pi*180) '°']; end
CombineImages(int,numel(focposv),numel(theta),'scale',s)

