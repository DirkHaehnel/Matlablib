% program simulating SOFI Fourier 

clear all
close all

NA = 1.45;
n0 = 1.51;
n = 1.46;
n1 = 1.46;
d0 = [];
d = 0;
d1 = [];
lamem = 0.625;
mag = 16/0.01;
pixel = 16;
pic = 1;
be_res = [];
al_res = [];
focusv = 0;
nn = 30;
mdf = MDFWideFieldMicroscope(NA,n0,n,n1,d0,d,d1,lamem,mag,focusv,pixel,pic,be_res,al_res,nn);
mdf = mdf/max(mdf(:));

Nmol = 100;
mm = 200+2*nn;
tst = zeros(mm);
for j=1:Nmol
    pos(j,:) = [ceil(mm*rand),ceil(mm*rand)];
    sta(j) = rand>=0.1;
end

kmax = 1000;
maxphot = 10;
back = 10;
res = zeros(mm,mm,kmax);
for k=1:kmax
    tst = zeros(mm,mm);
    for j=1:Nmol 
        if rand<0.01
            sta(j) = round(1-sta(j));
        end
        if sta(j)==1
            tst(pos(j,1),pos(j,2)) = 1;
        end
    end
    res(:,:,k) = poissrnd(back + maxphot*mConv2(tst,mdf));
    %     mim(res(:,:,k));
    %     caxis([0 5])
    %     eval(['print -dpng PalmSIM' mint2str(k,3)])
    k
end
res = res(nn+1:end-nn,nn+1:end-nn,:);

% ff = res;
% mm = max(res(:));
% for k=1:kmax
%     ff(:,:,k) = sum(exp(exp(i*(k-1)/kmax*2*pi)*res/mm),3);
%     k
% end
% ff = fft(log(ff),[],3);

meanim = mean(res,3);
delta = res-repmat(meanim,[1 1 size(res,3)]);
res = res./repmat(meanim,[1 1 size(res,3)])-1;
for m = 1:10
    k = 1; 
    cc(:,:,m,1) = res(:,:,k).*res(:,:,k+m); 
    for k=2:kmax-m 
        cc(:,:,m,1)=cc(:,:,m,1) + res(:,:,k).*res(:,:,k+m); 
    end
    cc(:,:,m,1) = cc(:,:,m,1)/(kmax-m);

    k = 1;
    cc(:,:,m,2) = res(:,:,k).*res(:,:,k+m).*res(:,:,k+2*m); 
    for k=2:kmax-2*m 
        cc(:,:,m,2)=cc(:,:,m,2) + res(:,:,k).*res(:,:,k+m).*res(:,:,k+2*m); 
    end
    cc(:,:,m,2) = cc(:,:,m,2)/(kmax-2*m);
    
    k = 1; 
    cc(:,:,m,3) = res(:,:,k).*res(:,:,k+m).*res(:,:,k+2*m).*res(:,:,k+3*m); 
    for k=2:kmax-3*m 
        cc(:,:,m,3) = cc(:,:,m,3) + res(:,:,k).*res(:,:,k+m).*res(:,:,k+2*m).*res(:,:,k+3*m); 
    end
    cc(:,:,m,3) = cc(:,:,m,3)/(kmax-3*m);

    k = 1;
    cc(:,:,m,4) = res(:,:,k).*res(:,:,k+m).*res(:,:,k+2*m).*res(:,:,k+3*m).*res(:,:,k+4*m); 
    for k=2:kmax-4*m 
        cc(:,:,m,4) = cc(:,:,m,4) + res(:,:,k).*res(:,:,k+m).*res(:,:,k+2*m).*res(:,:,k+3*m).*res(:,:,k+4*m); 
    end
    cc(:,:,m,4) = cc(:,:,m,4)/(kmax-4*m);
    
    k = 1; 
    cc(:,:,m,5) = res(:,:,k).*res(:,:,k+m).*res(:,:,k+2*m).*res(:,:,k+3*m).*res(:,:,k+4*m).*res(:,:,k+5*m); 
    for k=2:kmax-5*m 
        cc(:,:,m,5) = cc(:,:,m,5) + res(:,:,k).*res(:,:,k+m).*res(:,:,k+2*m).*res(:,:,k+3*m).*res(:,:,k+4*m).*res(:,:,k+5*m); 
    end
    cc(:,:,m,5) = cc(:,:,m,5)/(kmax-5*m);    
    
    %cc(:,:,m,3) = cc(:,:,m,3)-3*c(:,:,m,1);
    
%     -10 C[2] C[3] + DR[5]
%     -15 C[2]^3 - 10 C[3]^2 - 15 C[2] C[4] + DR[6]
%     -105 C[2]^2 C[3] - 35 C[3] C[4] - 21 C[2] C[5] + DR[7]
%     -105 C[2]^4 - 210 C[2]^2 C[4] - 35 C[4]^2 - 56 C[3] C[5] - 28 C[2] (10 C[3]^2 + C[6]) + DR[8]    
end
im4 = abs(cc(:,:,1,3)-cc(:,:,1,1).^2-cc(:,:,2,1).^2-cc(:,:,3,1).^2);
close; mim(cat(4,cat(3,meanim,cc(:,:,1,1)),cat(3,abs(cc(:,:,1,2)),im4)));

clear tst
tst(:,:,1) = mean(delta.^4-6*delta.^3+11*delta.^2,3)-3*mean(delta.^2,3).^2-6*meanim;
for k=2:4
    for j=1:size(delta,3)/2
        delta(:,:,j) = sum(delta(:,:,2*j-1:2*j),3);
    end
    delta = delta(:,:,1:end/2);
    tst(:,:,k) = mean(delta.^4-6*delta.^3+11*delta.^2,3)-3*mean(delta.^2,3).^2-6*meanim;
end

figure
m2 = (mm-2*nn-1)/2;
mdf = MDFWideFieldMicroscope(NA,n0,n,n1,d0,d,d1,lamem,mag,focusv,pixel,pic,be_res,al_res,m2);
otm2 = abs(fftshift(fft2(mdf.^2)));
otm4 = abs(fftshift(fft2(mdf.^4)));
otm2 = otm2/max(otm2(:));
otm4 = otm4/max(otm4(:));
otm2a = interp2((-2*m2:2:2*m2)',-2*m2:2:2*m2,otm1,(-m2:m2)',-m2:m2,'cubic');
otm2a = otm2a/max(otm2a(:));
otm4a = interp2((-4*m2:4:4*m2)',-4*m2:4:4*m2,otm1,(-m2:m2)',-m2:m2,'cubic');
otm4a = otm4a/max(otm4a(:));

weight2 = otm2a./(0.01+otm2).*(otm2>0.01);
weight4 = otm4a./(0.01+otm4).*(otm4>0.01);
tst2 = abs(ifft2(ifftshift(weight2.*fftshift(fft2(cc(:,:,1,1))))));
tst4 = abs(ifft2(ifftshift(weight4.*fftshift(fft2(im4)))));
clf
mim(cat(4,cat(3,cc(:,:,1,1),abs(cc(:,:,1,3)-cc(:,:,1,1).^2-cc(:,:,2,1).^2-cc(:,:,3,1).^2)),cat(3,tst2,tst4)));

bla=interp2((-m2:m2)',-m2:m2,mdf,(-4*m2:4:4*m2)',-4*m2:4:4*m2,'cubic',0);
mim(cat(3,im4,abs(mConv2(im4,bla-0.666*mdf.^4))))