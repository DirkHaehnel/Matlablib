% calculates a simplified estimate of the DIC PSF

lamem = 0.67;
lamex = 0.64;
n = 1.33;
NA = 1.2;
w0 = 0.4;
av = 50/60;
delta = 0.4;

kem=2*pi/lamem; 
z=(0:0.1:2);

d = delta/2;
w0 = LaserBeamFit(w0,lamex/n,z);
[x,y] = meshgrid(-1:1/100:1,1/200:1/100:1);
rr = sqrt(x.^2+y.^2);
pp = angle(x+i*y);
mm = ones(size(x));
for j=1:size(z,2)
    wm = 3*max(w0(j));
    res = D3CEF(wm*rr*kem,mm*z(j)*kem,av*kem,NA,n).*exp(-2*((wm*x-d).^2+(wm*y).^2)/w0(j).^2);    
    nrm = sum(res(:));
    mx(j) = wm*sum(sum(x.*res))/nrm;
    wx(j) = 2*sqrt(sum(sum((wm*x-mx(j)).^2.*res))/nrm);
    wy(j) = 2*sqrt(sum(sum((wm*y).^2.*res))/nrm);
    amp(j) = D3CEF(0,z(j)*kem,av*kem,NA,n);    
    amp1(j) = nrm*w0(j).^2;    
end   

