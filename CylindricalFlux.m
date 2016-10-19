function [flux, fm, fn] = CylindricalFlux(w,mv,cjm,cjn,chm,chn,n,r)

% CylindricalFlux calculates the energy flux through a cylindrical surface with
% radius r if given the cylindrical coefficients cjm, cjn, chm,chn 
% of the electric field and the refractive index n

if 1
    bj = inline('besselj(n,x,1)','n','x');
    bh = inline('besselh(n,1,x,1)','n','x');
    bjx = inline('0.5*(besselj(n-1,x,1) - besselj(n+1,x,1))','n','x');
    bhx = inline('0.5*(besselh(n-1,1,x,1) - besselh(n+1,1,x,1))','n','x');
else
    bj = inline('besselj(n,x)','n','x');
    bh = inline('besselh(n,x)','n','x');
    bjx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x,1))','n','x');
    bhx = inline('0.5*(besselh(n-1,x) - besselh(n+1,x))','n','x');
end

w = w(:);
q = sqrt(n^2-w.^2);
flux = 0;
for j=1:length(mv)
    m = mv(j);
    fm = imag((chm(:,j).*bhx(m,q*r)+cjm(:,j).*bjx(m,q*r)).*q.*conj(q.^2.*(chm(:,j).*bh(m,q*r)+cjm(:,j).*bj(m,q*r))));
    flux = flux + sum(fm(isfinite(fm)));
    fn = imag(q.^2.*(chn(:,j).*bh(m,q*r)+cjn(:,j).*bj(m,q*r)).*conj(q*n^2.*(chn(:,j).*bhx(m,q*r)+cjn(:,j).*bjx(m,q*r))));
    flux = flux - sum(fn(isfinite(fn)));
end 
flux = flux*pi*r*mean(diff(w));

