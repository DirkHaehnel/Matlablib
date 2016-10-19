% program for the calculation of the change in extinction

load metals
warning off

bj = inline('besselj(n+0.5,x).*sqrt(pi./x/2)','n','x');
bh = inline('besselh(n+0.5,x).*sqrt(pi./x/2)','n','x');
bjx = inline('sqrt(pi/8./x).*(besselj(n+0.5,x)./x + besselj(n-0.5,x) - besselj(n+1.5,x))','n','x');
bhx = inline('sqrt(pi/8./x).*(besselh(n+0.5,x)./x + besselh(n-0.5,x) - besselh(n+1.5,x))','n','x');

kout = 1.;
%kout = 2.5;
d = 50; 
lambda = 560;
N = 30; n = 1:N;
rm = d + (0.5:100)/100*25;

%kin = 1; 
kin = silver(wavelength==lambda); 
rv = d/lambda*2*pi;
nv = [kin kout];
resx = []; resz = [];      
for k = 1:length(rm)
    r = rm(k)/lambda*2*pi;
    
    % coefficients of the dipole fields, as required by Spherical
    % x-orientation
    cx = [-kout*sqrt(pi*(2*n+1)./n./(n+1)).*bh(n,kout*r); ...
            -1i*kout*sqrt(pi*(2*n+1)./n./(n+1)).*bhx(n,kout*r); zeros(2,N)];
    
    % z-orientation
    cz = [zeros(1,N); 1i*sqrt(4*pi*(2*n+1)).*bh(n,kout*r)/r; zeros(2,N)];        
    
    coefx = Spherical(n,nv,rv,cx);
    coefz = Spherical(n,nv,rv,cz);
    
    % emission intensity normalized to polymer
    cxh = [-kout*sqrt(pi*(2*n+1)./n./(n+1)).*bj(n,kout*r); ...
            -1i*kout*sqrt(pi*(2*n+1)./n./(n+1)).*bjx(n,kout*r)];
    resx(k,2) = 2*SphericalFlux(zeros(1,N),zeros(1,N),coefx(1,:)+cxh(1,:),coefx(2,:)+cxh(2,:),kout,inf); 
    resx(k,1) = resx(k,2) - 2*SphericalFlux(cx(1,:),cx(2,:),coefx(1,:),coefx(2,:),kout,rv);
    czh = 1i*sqrt(4*pi*(2*n+1)).*bj(n,kout*r)/r;        
    resz(k,2) = SphericalFlux(zeros(1,N),zeros(1,N),zeros(1,N),czh+coefz(2,:),kout,inf);
    resz(k,1) = resz(k,2) - SphericalFlux(cz(1,:),cz(2,:),coefz(1,:),coefz(2,:),kout,rv);                    
end    

c0 = 2*pi*1i.^(n+1).*sqrt((2*n+1)/4/pi./n./(n+1));
nn = (n'.*(n+1)')*ones(1,length(rm));
r = rm/lambda*2*pi;
c = [c0; c0; zeros(2,N)];
coef = Spherical(n,nv,rv,c);
%dm = sqrt(imag(conj(coef(1,:).*bh(n,kout*rv) + c(1,:).*bj(n,kout*rv)).*(coef(1,:).*bhx(n,kout*rv) + c(1,:).*bjx(n,kout*rv))))*kout*rv;
%dn = sqrt(imag(conj(coef(2,:).*bh(n,kout*rv) + c(2,:).*bj(n,kout*rv)).*(coef(2,:).*bhx(n,kout*rv) + c(2,:).*bjx(n,kout*rv))))*kout*rv;
for k=1:N
    cm(k,:) = coef(1,k)*bh(k,kout*r) + c(1,k)*bj(k,kout*r);
    cn(k,:) = coef(2,k)*bhx(k,kout*r) + c(2,k)*bjx(k,kout*r);
    cr(k,:) = coef(2,k)*bh(k,kout*r)./r/kout + c(2,k)*bj(k,kout*r)./r/kout;
    dm(k,:) = sqrt(imag(conj(coef(1,k)*bh(k,kout*rv) + c(1,k)*bj(k,kout*rv)).*...
        (coef(1,k)*bhx(k,kout*rv) + c(1,k)*bjx(k,kout*rv))))*kout*rv*bh(k,kout*r);
    dn(k,:) = sqrt(imag(conj(coef(2,k)*bh(k,kout*rv) + c(2,k)*bj(k,kout*rv)).*...
        (coef(2,k)*bhx(k,kout*rv) + c(2,k)*bjx(k,kout*rv))))*kout*rv*bhx(k,kout*r);
    dr(k,:) = sqrt(imag(conj(coef(2,k)*bh(k,kout*rv) + c(2,k)*bj(k,kout*rv)).*...
        (coef(2,k)*bhx(k,kout*rv) + c(2,k)*bjx(k,kout*rv))))*kout*rv*bh(k,kout*r)./r/kout;
end
intx = sum(nn.*(abs(cm).^2 + abs(cn).^2))/4/pi*kout;
intz = sum(abs(nn.*cr).^2)/2/pi*kout;
flx = sum(nn.*(abs(dm).^2 + abs(dn).^2))/4/pi*kout;
flz = sum(abs(nn.*dr).^2)/2/pi*kout;

semilogy(rm,resx(:,1)'./(intx+flx),rm,resz(:,1)'./(intz+flz))