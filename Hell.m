function z = Hell(y,t,tau,kisc,kph)

% function Hell calculates a series of convoluiton kernels for the
% intensity profile y
% (c) 2004 Joerg Enderlein

if size(y,1)==1 || size(y,2)==1
    y = y(:);
    f = y./(1+tau*y);
    f0 = kph*f./(kph+tau*kisc*f);
    fac = tau*kisc*f.^2./(kph+tau*kisc*f);
    if length(t)==2
        z = f0*diff(t) + fac.*(exp(-(kph+tau*kisc*f)*t(1))-exp(-(kph+tau*kisc*f)*t(2)))./(kph+tau*kisc*f);
    else
        for j= 1:length(t)
            z(:,j) = f0 + fac.*exp(-(kph+tau*kisc*f)*t(j));
        end
    end
else
    f = y./(1+tau*y);
    f0 = kph*f./(kph+tau*kisc*f);
    fac = tau*kisc*f.^2./(kph+tau*kisc*f);
    if length(t)==2
        z = f0*diff(t) + fac.*(exp(-(kph+tau*kisc*f)*t(1))-exp(-(kph+tau*kisc*f)*t(2)))./(kph+tau*kisc*f);
    else
        for j= 1:length(t)
            z(:,:,j) = f0 + fac.*exp(-(kph+tau*kisc*f)*t(j));
        end
    end
end

return

im = sum(double(imread('D:\Joerg\Doc\Hell\LeeuwenhoekBig.jpg')),3);
im0 = mconv2(im,z1);

exc0 = 1e7;
tau = 1e-9;
kisc = 1e8;
kph = 1e5;
tv = [0.1 0.2 0.4 0.8 1.6 3.2 6.4]*1e-6;
tauv = [1 2 4 8]*1e-6;
t = [0 tv(1)];
clear imm
for j=2:length(tv)
    t = [tv(j-1) tv(j)];
    imm(:,:,j-1) = mconv2(im,Hell(exc0*z1,t,tau,kisc,kph));
    MM(j-1,:) = [diff(t) -diff(exp(-t'*(1./tauv))).*tauv];
end
imm = reshape(imm,size(imm,1)*size(imm,2),size(imm,3));
tmp = inv(MM'*MM)*MM'*imm';
im1 = reshape(tmp(2,:).^2./sqrt(sum(tmp(2:end,:).^2)),size(im,1),size(im,2));

for j=jj:size(imm,1)
    bla(:,j) = lsqnonneg(MM,imm(j,:)');
end

clear imm
for j=2:length(tv)
    t = [tv(j-1) tv(j)];
    imm(:,:,j-1) = mconv2(im,Hell(exc0*z2,t,tau,kisc,kph));
end
imm = reshape(imm,size(imm,1)*size(imm,2),size(imm,3));
tmp = inv(MM'*MM)*MM'*imm';
im2 = reshape(tmp(2,:).^2./sqrt(sum(tmp(2:end,:).^2)),size(im,1),size(im,2));

for j=jj:size(imm,1)
    bla(:,j) = lsqnonneg(MM,imm(j,:)');
end

tst=im;
tst(90:130,end-500:end-496)=0;
tst(90:130,end-104:end-100)=0;
tst(100:120,end-500:end-100)=0;

subplot(2,2,1); mim(tst); subplot(2,2,2); mim(im0); subplot(2,2,3); mim(sqrt(im1)); subplot(2,2,4); mim(sqrt(im2));


return

NA = 1.4;
n = 1.51;
over = [NA*3e3 inf 5e3];
lamex = 0.514;
[fx0, fx2, fz, rho, z, exc, excx, excz] = NewExc([-0.00125 2], 0, NA, n, n, n, [], 0, [], lamex, over, 0, [], lamex/0.0025);
exc = exc(:,1)';
rho = rho(:,1)';
clear fx0 fx2 fz z excx ecxz

y = -1.5:0.0025:1.5;
[x,y] = meshgrid(y,y);
x = sqrt(x.^2+y.^2);
z2 = cos(asin(NA/n)).^2*besselj(1,2*pi*NA/lamex*x).^2+(NA/n).^2*besselj(0,2*pi*NA/lamex*x).^2;
z2 = z2.*disk((size(z2,1)-1)/2);
z2 = z2/max(z2(:));
z1 = interp1(rho,exc,x,'cubic');
clear x y

% for j=1:size(im,1) for k=1:size(im,2) im00(j,k) = sum(sum(im.*z2([j:-1:2 1:end-j+1],[k:-1:2 1:end-k+1]))); end; end

return

% PSF calculation

y = -1.:0.0025:1.;
im = zeros(length(y));
im(:,(size(im,2)-1)/2) = 1;
im0 = mconv2(im,z1);
clear imm
exc0 = 1e7;
tau = 1e-9;
kisc = 1e8;
kph = 1e5;
tv = [0.1 0.2 0.4 0.8 1.6 3.2 6.4]*1e-6;
tauv = [1 2 4 8]*1e-6;
t = [0 tv(1)];

for j=2:length(tv)
    t = [tv(j-1) tv(j)];
    imm(:,:,j-1) = mconv2(im,Hell(exc0*z1,t,tau,kisc,kph));
    MM(j-1,:) = [diff(t) -diff(exp(-t'*(1./tauv))).*tauv];
end
for j=1:size(im,2) tst1(:,j)=lsqnonneg(MM,squeeze(imm(k,j,:))); end

clear imm
for j=2:length(tv)
    t = [tv(j-1) tv(j)];
    imm(:,:,j-1) = mconv2(im,Hell(exc0*z2,t,tau,kisc,kph));
    MM(j-1,:) = [diff(t) -diff(exp(-t'*(1./tauv))).*tauv];
end
for j=1:size(im,2) tst2(:,j)=lsqnonneg(MM,squeeze(imm(1,j,:))); end

plot(y,im0((end-1)/2,:)/max(im0(:)),y,tst1(2,:)/max(tst1(2,:)),y,tst2(2,:)/max(tst2(2,:)),[-1 1],[1 1]/exp(2),':')
axis([-0.75 0.75 0 1])
xlabel('\itx\rm [\mum]');
ylabel('rel. amplitude');

sqrt([sum(y.^2.*im0((end-1)/2,:))/sum(im0((end-1)/2,:)),sum(y.^2.*tst1(2,:))/sum(tst1(2,:)),sum(y.^2.*tst2(2,:))/sum(tst2(2,:))])
