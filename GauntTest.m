% GauntTest

if 0
    
    x = -1:1e-3:1;

    n1 = 80;
    n2 = 20;
    m1 = 8;
    m2 = 12;

    tst = legendre(n1,x,'norm');
    f0 = tst(m1+1,:);
    tst = legendre(n2,x,'norm');
    f0 = f0.*tst(m2+1,:);

    fac = gaunt(n1+n2,n1,n2,m1+m2,m1,m2,'recursive');
    f1 = 0*x;
    for j=max([m1+m2,abs(n1-n2)]):(n1+n2)
        tst = legendre(j,x,'norm');
        f1 = f1 + fac(j+1)*tst(m1+m2+1,:);
    end

    plot(x,f0,x,f1*2.5)

end

if 1

    theta = 0:pi/360:pi;
    phi = 0;
    rad = (0.05:0.1:10)';

    n0 = 10;
    m0 = 5;

    tst = legendre(n0,cos(theta),'norm');

    f0 = (sqrt(pi/2./rad).*besselj(n0+0.5,rad))*tst(m0+1,:);

    % surf(rad*sin(theta),rad*cos(theta),f0); shading interp; axis image

    dd = 5;
    tt = pi/6;

    nmax = 20;

    bes = sqrt(pi/2./dd).*besselj((0:n0+nmax)+0.5,dd);
    for n=1:nmax+n0
        M{n} = legendre(n,cos(tt),'norm');
    end
    f1 = 0*f0;
    for n=1:nmax
        leg = legendre(n,cos(theta),'norm');
        for m=-n:n
            g = gaunt(n0+n,n0,n,m0-m,m0,-m,'recursive');
            p = max([abs(m),abs(n0-n)]):(n0+n);
            clear tmp;
            for k=1:length(p)
                if abs(m0-m)<=p(k)
                    tmp(k) = M{p(k)}(abs(m0-m)+1);
                else
                    tmp(k) = 0;
                end
            end
            coef = (-1)^m*(2*n+1)*i^(n-n0)*sum((-i).^p.*(2*p+1).*g(p+1).*bes(p+1).*tmp);
            f1 = f1 + coef*(sqrt(pi/2./rad).*besselj(n+0.5,rad))*leg(abs(m)+1,:);
        end
    end
    subplot(121); pcolor(rad*sin(theta),rad*cos(theta),f0);shading interp; axis image; 
    subplot(122); pcolor(rad*sin(theta),rad*cos(theta),f1);shading interp; axis image
end


