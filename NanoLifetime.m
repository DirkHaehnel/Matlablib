% program for the calculation of fluorescence around a metallized sphere

load metals
warning off

bj = inline('besselj(n+0.5,x).*sqrt(pi./x/2)','n','x');
bh = inline('besselh(n+0.5,x).*sqrt(pi./x/2)','n','x');
bjx = inline('sqrt(pi/8./x)*(besselj(n+0.5,x)./x + besselj(n-0.5,x) - besselj(n+1.5,x))','n','x');
bhx = inline('sqrt(pi/8./x)*(besselh(n+0.5,x)./x + besselh(n-0.5,x) - besselh(n+1.5,x))','n','x');

kbead = 1.5;
kwater = 1.333;
d = 5; 
rov = (22.5:0.5:100); 
lamem = [560 670];
lamex = [514 635];
N = 40; n = 1:N;
rm = (0.5:130)/130*100;
    
qux = zeros(length(rm),length(rov),length(lamem)); 
quz = qux;
emx = qux; 
emz = qux;
ltx = qux;
ltz = qux;

load silverlifetime 
for s=2:length(lamem)
    kem = silver(wavelength==lamem(s));
    for j = 1:length(rov)
        disp(['em = ' num2str(lamem(s)) ' & j = ' num2str(j)]);
        if rov(j)<=d
            rv = rov(j)/lamem(s)*2*pi;
            nv = [kem kwater];
        else            
            rv = [rov(j)-d rov(j)]/lamem(s)*2*pi;
            nv = [kbead kem kwater];
        end
        resx = []; resz = [];      
        for k = 1:length(rm)
            r = rm(k)/lamem(s)*2*pi;
            
            % coefficients of the dipole fields, as required by Spherical
            % x-orientation
            if r<rv(1)
                cx = [zeros(4*length(rv)-2,N); -nv(1)*sqrt(pi*(2*n+1)./n./(n+1)).*bj(n,nv(1)*r); ...
                        -i*nv(1)*sqrt(pi*(2*n+1)./n./(n+1)).*bjx(n,nv(1)*r)];
            else                
                cx = [-nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bh(n,nv(end)*r); ...
                        -i*nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bhx(n,nv(end)*r); zeros(4*length(rv)-2,N)];
            end 

            % z-orientation
            if r<rv(1)    
                cz = [zeros(4*length(rv)-1,N); i*sqrt(4*pi*(2*n+1)).*bj(n,nv(1)*r)/r];
            else
                cz = [zeros(1,N); i*sqrt(4*pi*(2*n+1)).*bh(n,nv(end)*r)/r; zeros(4*length(rv)-2,N)];        
            end

            if rm(k)<rov(j)-d | r>rv(end)
                coefx = Spherical(n,nv,rv,cx);
                coefz = Spherical(n,nv,rv,cz);
                
                % emission intensity normalized to polymer
                if r<rv(1)
                    resx(k,1) = 2*SphericalFlux(coefx(end-1,:),coefx(end,:),cx(end-1,:),cx(end,:),nv(1),rv(1));
                    resx(k,2) = 2*SphericalFlux(zeros(1,N),zeros(1,N),coefx(1,:),coefx(2,:),nv(end),rv(end));
                    resz(k,1) = SphericalFlux(coefz(end-1,:),coefz(end,:),cz(end-1,:),cz(end,:),nv(1),rv(1));
                    resz(k,2) = SphericalFlux(zeros(1,N),zeros(1,N),coefz(1,:),coefz(2,:),nv(end),rv(end));
                else
                    cxh = [-nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bj(n,nv(end)*r); ...
                        -i*nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bjx(n,nv(end)*r)];
                    resx(k,2) = 2*SphericalFlux(zeros(1,N),zeros(1,N),coefx(1,:)+cxh(1,:),coefx(2,:)+cxh(2,:),nv(end),inf);
                    resx(k,1) = resx(k,2) - 2*SphericalFlux(cx(1,:),cx(2,:),coefx(1,:),coefx(2,:),nv(end),rv(end));
                    czh = i*sqrt(4*pi*(2*n+1)).*bj(n,nv(end)*r)/r;        
                    resz(k,2) = SphericalFlux(zeros(1,N),zeros(1,N),zeros(1,N),czh+coefz(2,:),nv(end),inf);
                    resz(k,1) = resz(k,2) - SphericalFlux(cz(1,:),cz(2,:),coefz(1,:),coefz(2,:),nv(end),rv(end));                    
                end        
            else
                resx(k,1:2) = inf;
                resz(k,1:2) = inf;
            end
        end    
        qux(:,j,s) = resx(:,2)./resx(:,1); 
        emx(:,j,s) = resx(:,2)/(kwater/3);
        ltx(:,j,s) = kwater/3./resx(:,1);
        quz(:,j,s) = resz(:,2)./resz(:,1); 
        emz(:,j,s) = resz(:,2)/(kwater/3);
        ltz(:,j,s) = kwater/3./resz(:,1);
    end
end
save silverlifetime qux quz emx emz ltx ltz rov rm lamex lamem N kbead kwater d
 
c0 = 2*pi*i.^(n+1).*sqrt((2*n+1)/4/pi./n./(n+1));
nn = (n'.*(n+1)')*ones(1,length(rm));
for s=1:length(lamex)
    kex = silver(wavelength==lamex(s));
    nv = [kbead kex kwater];
    r = rm/lamex(s)*2*pi;
    for j = 1:length(rov)
        disp(['ex = ' num2str(lamex(s)) ' & j = ' num2str(j)]);
        if rov(j)<=d
            rv = rov(j)/lamex(s)*2*pi;
            c = [c0; c0; zeros(4*length(rv)-2,N)];
            coef = Spherical(n,nv(2:end),rv,c);
        else            
            rv = [rov(j)-d rov(j)]/lamex(s)*2*pi;
            c = [c0; c0; zeros(4*length(rv)-2,N)];
            coef = Spherical(n,nv,rv,c);
        end
        ind1 = rm<rov(j)-d; ind2 = r>rv(end); ind3 = ~(ind1 | ind2);        
        for k=1:N
            cm(k,ind1) = coef(end-1,k)*bj(k,nv(1)*r(ind1));
            cm(k,ind3) = 0;
            cm(k,ind2) = coef(1,k)*bh(k,nv(end)*r(ind2)) + c(1,k)*bj(k,nv(end)*r(ind2));
            cn(k,ind1) = coef(end,k)*bjx(k,nv(1)*r(ind1));
            cn(k,ind3) = 0;
            cn(k,ind2) = coef(2,k)*bhx(k,nv(end)*r(ind2)) + c(1,k)*bjx(k,nv(end)*r(ind2));
            cr(k,ind1) = coef(end,k)*bj(k,nv(1)*r(ind1))./r(ind1)/nv(1);
            cr(k,ind3) = 0;
            cr(k,ind2) = coef(2,k)*bh(k,nv(end)*r(ind2))./r(ind2)/nv(end) + c(1,k)*bj(k,nv(end)*r(ind2))./r(ind2)/nv(end);
        end
        intx(j,:,s) = sum(nn.*(abs(cm).^2 + abs(cn).^2))/4/pi;
        intz(j,:,s) = sum(abs(nn.*cr).^2)/2/pi;
    end
end
save silverlifetime intx intz -append


qux = zeros(length(rm),length(rov),length(lamem)); 
quz = qux;
emx = qux; 
emz = qux;
ltx = qux;
ltz = qux;

for s=1:length(lamem)
    kem = gold(wavelength==lamem(s));
    for j = 1:length(rov)
        disp(['em = ' num2str(lamem(s)) ' & j = ' num2str(j)]);
        if rov(j)<=d
            rv = rov(j)/lamem(s)*2*pi;
            nv = [kem kwater];
        else            
            rv = [rov(j)-d rov(j)]/lamem(s)*2*pi;
            nv = [kbead kem kwater];
        end
        resx = []; resz = [];      
        for k = 1:length(rm)
            r = rm(k)/lamem(s)*2*pi;
            
            % coefficients of the dipole fields, as required by Spherical
            % x-orientation
            if r<rv(1)
                cx = [zeros(4*length(rv)-2,N); -nv(1)*sqrt(pi*(2*n+1)./n./(n+1)).*bj(n,nv(1)*r); ...
                        -i*nv(1)*sqrt(pi*(2*n+1)./n./(n+1)).*bjx(n,nv(1)*r)];
            else                
                cx = [-nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bh(n,nv(end)*r); ...
                        -i*nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bhx(n,nv(end)*r); zeros(4*length(rv)-2,N)];
            end 

            % z-orientation
            if r<rv(1)    
                cz = [zeros(4*length(rv)-1,N); i*sqrt(4*pi*(2*n+1)).*bj(n,nv(1)*r)/r];
            else
                cz = [zeros(1,N); i*sqrt(4*pi*(2*n+1)).*bh(n,nv(end)*r)/r; zeros(4*length(rv)-2,N)];        
            end

            if rm(k)<rov(j)-d | r>rv(end)
                coefx = Spherical(n,nv,rv,cx);
                coefz = Spherical(n,nv,rv,cz);
                
                % emission intensity normalized to polymer
                if r<rv(1)
                    resx(k,1) = 2*SphericalFlux(coefx(end-1,:),coefx(end,:),cx(end-1,:),cx(end,:),nv(1),rv(1));
                    resx(k,2) = 2*SphericalFlux(zeros(1,N),zeros(1,N),coefx(1,:),coefx(2,:),nv(end),rv(end));
                    resz(k,1) = SphericalFlux(coefz(end-1,:),coefz(end,:),cz(end-1,:),cz(end,:),nv(1),rv(1));
                    resz(k,2) = SphericalFlux(zeros(1,N),zeros(1,N),coefz(1,:),coefz(2,:),nv(end),rv(end));
                else
                    cxh = [-nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bj(n,nv(end)*r); ...
                        -i*nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bjx(n,nv(end)*r)];
                    resx(k,2) = 2*SphericalFlux(zeros(1,N),zeros(1,N),coefx(1,:)+cxh(1,:),coefx(2,:)+cxh(2,:),nv(end),inf);
                    resx(k,1) = resx(k,2) - 2*SphericalFlux(cx(1,:),cx(2,:),coefx(1,:),coefx(2,:),nv(end),rv(end));
                    czh = i*sqrt(4*pi*(2*n+1)).*bj(n,nv(end)*r)/r;        
                    resz(k,2) = SphericalFlux(zeros(1,N),zeros(1,N),zeros(1,N),czh+coefz(2,:),nv(end),inf);
                    resz(k,1) = resz(k,2) - SphericalFlux(cz(1,:),cz(2,:),coefz(1,:),coefz(2,:),nv(end),rv(end));                    
                end        
            else
                resx(k,1:2) = inf;
                resz(k,1:2) = inf;
            end
        end    
        qux(:,j,s) = resx(:,2)./resx(:,1); 
        emx(:,j,s) = resx(:,2)/(kwater/3);
        ltx(:,j,s) = kwater/3./resx(:,1);
        quz(:,j,s) = resz(:,2)./resz(:,1); 
        emz(:,j,s) = resz(:,2)/(kwater/3);
        ltz(:,j,s) = kwater/3./resz(:,1);
    end
end
save goldlifetime qux quz emx emz ltx ltz rov rm lamex lamem N kbead kwater d
 
c0 = 2*pi*i.^(n+1).*sqrt((2*n+1)/4/pi./n./(n+1));
nn = (n'.*(n+1)')*ones(1,length(rm));
for s=1:length(lamex)
    kex = gold(wavelength==lamex(s));
    nv = [kbead kex kwater];
    r = rm/lamex(s)*2*pi;
    for j = 1:length(rov)
        disp(['ex = ' num2str(lamex(s)) ' & j = ' num2str(j)]);
        if rov(j)<=d
            rv = rov(j)/lamex(s)*2*pi;
            c = [c0; c0; zeros(4*length(rv)-2,N)];
            coef = Spherical(n,nv(2:end),rv,c);
        else            
            rv = [rov(j)-d rov(j)]/lamex(s)*2*pi;
            c = [c0; c0; zeros(4*length(rv)-2,N)];
            coef = Spherical(n,nv,rv,c);
        end
        ind1 = rm<rov(j)-d; ind2 = r>rv(end); ind3 = ~(ind1 | ind2);        
        for k=1:N
            cm(k,ind1) = coef(end-1,k)*bj(k,nv(1)*r(ind1));
            cm(k,ind3) = 0;
            cm(k,ind2) = coef(1,k)*bh(k,nv(end)*r(ind2)) + c(1,k)*bj(k,nv(end)*r(ind2));
            cn(k,ind1) = coef(end,k)*bjx(k,nv(1)*r(ind1));
            cn(k,ind3) = 0;
            cn(k,ind2) = coef(2,k)*bhx(k,nv(end)*r(ind2)) + c(1,k)*bjx(k,nv(end)*r(ind2));
            cr(k,ind1) = coef(end,k)*bj(k,nv(1)*r(ind1))./r(ind1)/nv(1);
            cr(k,ind3) = 0;
            cr(k,ind2) = coef(2,k)*bh(k,nv(end)*r(ind2))./r(ind2)/nv(end) + c(1,k)*bj(k,nv(end)*r(ind2))./r(ind2)/nv(end);
        end
        intx(j,:,s) = sum(nn.*(abs(cm).^2 + abs(cn).^2))/4/pi;
        intz(j,:,s) = sum(abs(nn.*cr).^2)/2/pi;
    end
end
save goldlifetime intx intz -append


qux = zeros(length(rm),length(rov),length(lamem)); 
quz = qux;
emx = qux; 
emz = qux;
ltx = qux;
ltz = qux;

for s=1:length(lamem)
    kem = platinum(wavelength==lamem(s));
    for j = 1:length(rov)
        disp(['em = ' num2str(lamem(s)) ' & j = ' num2str(j)]);
        if rov(j)<=d
            rv = rov(j)/lamem(s)*2*pi;
            nv = [kem kwater];
        else            
            rv = [rov(j)-d rov(j)]/lamem(s)*2*pi;
            nv = [kbead kem kwater];
        end
        resx = []; resz = [];      
        for k = 1:length(rm)
            r = rm(k)/lamem(s)*2*pi;
            
            % coefficients of the dipole fields, as required by Spherical
            % x-orientation
            if r<rv(1)
                cx = [zeros(4*length(rv)-2,N); -nv(1)*sqrt(pi*(2*n+1)./n./(n+1)).*bj(n,nv(1)*r); ...
                        -i*nv(1)*sqrt(pi*(2*n+1)./n./(n+1)).*bjx(n,nv(1)*r)];
            else                
                cx = [-nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bh(n,nv(end)*r); ...
                        -i*nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bhx(n,nv(end)*r); zeros(4*length(rv)-2,N)];
            end 

            % z-orientation
            if r<rv(1)    
                cz = [zeros(4*length(rv)-1,N); i*sqrt(4*pi*(2*n+1)).*bj(n,nv(1)*r)/r];
            else
                cz = [zeros(1,N); i*sqrt(4*pi*(2*n+1)).*bh(n,nv(end)*r)/r; zeros(4*length(rv)-2,N)];        
            end

            if rm(k)<rov(j)-d | r>rv(end)
                coefx = Spherical(n,nv,rv,cx);
                coefz = Spherical(n,nv,rv,cz);
                
                % emission intensity normalized to polymer
                if r<rv(1)
                    resx(k,1) = 2*SphericalFlux(coefx(end-1,:),coefx(end,:),cx(end-1,:),cx(end,:),nv(1),rv(1));
                    resx(k,2) = 2*SphericalFlux(zeros(1,N),zeros(1,N),coefx(1,:),coefx(2,:),nv(end),rv(end));
                    resz(k,1) = SphericalFlux(coefz(end-1,:),coefz(end,:),cz(end-1,:),cz(end,:),nv(1),rv(1));
                    resz(k,2) = SphericalFlux(zeros(1,N),zeros(1,N),coefz(1,:),coefz(2,:),nv(end),rv(end));
                else
                    cxh = [-nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bj(n,nv(end)*r); ...
                        -i*nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bjx(n,nv(end)*r)];
                    resx(k,2) = 2*SphericalFlux(zeros(1,N),zeros(1,N),coefx(1,:)+cxh(1,:),coefx(2,:)+cxh(2,:),nv(end),inf);
                    resx(k,1) = resx(k,2) - 2*SphericalFlux(cx(1,:),cx(2,:),coefx(1,:),coefx(2,:),nv(end),rv(end));
                    czh = i*sqrt(4*pi*(2*n+1)).*bj(n,nv(end)*r)/r;        
                    resz(k,2) = SphericalFlux(zeros(1,N),zeros(1,N),zeros(1,N),czh+coefz(2,:),nv(end),inf);
                    resz(k,1) = resz(k,2) - SphericalFlux(cz(1,:),cz(2,:),coefz(1,:),coefz(2,:),nv(end),rv(end));                    
                end        
            else
                resx(k,1:2) = inf;
                resz(k,1:2) = inf;
            end
        end    
        qux(:,j,s) = resx(:,2)./resx(:,1); 
        emx(:,j,s) = resx(:,2)/(kwater/3);
        ltx(:,j,s) = kwater/3./resx(:,1);
        quz(:,j,s) = resz(:,2)./resz(:,1); 
        emz(:,j,s) = resz(:,2)/(kwater/3);
        ltz(:,j,s) = kwater/3./resz(:,1);
    end
end
save platinumlifetime qux quz emx emz ltx ltz rov rm lamex lamem N kbead kwater d
 
c0 = 2*pi*i.^(n+1).*sqrt((2*n+1)/4/pi./n./(n+1));
nn = (n'.*(n+1)')*ones(1,length(rm));
for s=1:length(lamex)
    kex = platinum(wavelength==lamex(s));
    nv = [kbead kex kwater];
    r = rm/lamex(s)*2*pi;
    for j = 1:length(rov)
        disp(['ex = ' num2str(lamex(s)) ' & j = ' num2str(j)]);
        if rov(j)<=d
            rv = rov(j)/lamex(s)*2*pi;
            c = [c0; c0; zeros(4*length(rv)-2,N)];
            coef = Spherical(n,nv(2:end),rv,c);
        else            
            rv = [rov(j)-d rov(j)]/lamex(s)*2*pi;
            c = [c0; c0; zeros(4*length(rv)-2,N)];
            coef = Spherical(n,nv,rv,c);
        end
        ind1 = rm<rov(j)-d; ind2 = r>rv(end); ind3 = ~(ind1 | ind2);        
        for k=1:N
            cm(k,ind1) = coef(end-1,k)*bj(k,nv(1)*r(ind1));
            cm(k,ind3) = 0;
            cm(k,ind2) = coef(1,k)*bh(k,nv(end)*r(ind2)) + c(1,k)*bj(k,nv(end)*r(ind2));
            cn(k,ind1) = coef(end,k)*bjx(k,nv(1)*r(ind1));
            cn(k,ind3) = 0;
            cn(k,ind2) = coef(2,k)*bhx(k,nv(end)*r(ind2)) + c(1,k)*bjx(k,nv(end)*r(ind2));
            cr(k,ind1) = coef(end,k)*bj(k,nv(1)*r(ind1))./r(ind1)/nv(1);
            cr(k,ind3) = 0;
            cr(k,ind2) = coef(2,k)*bh(k,nv(end)*r(ind2))./r(ind2)/nv(end) + c(1,k)*bj(k,nv(end)*r(ind2))./r(ind2)/nv(end);
        end
        intx(j,:,s) = sum(nn.*(abs(cm).^2 + abs(cn).^2))/4/pi;
        intz(j,:,s) = sum(abs(nn.*cr).^2)/2/pi;
    end
end
save platinumlifetime intx intz -append


qux = zeros(length(rm),length(rov),length(lamem)); 
quz = qux;
emx = qux; 
emz = qux;
ltx = qux;
ltz = qux;

for s=1:length(lamem)
    kem = copper(wavelength==lamem(s));
    for j = 1:length(rov)
        disp(['em = ' num2str(lamem(s)) ' & j = ' num2str(j)]);
        if rov(j)<=d
            rv = rov(j)/lamem(s)*2*pi;
            nv = [kem kwater];
        else            
            rv = [rov(j)-d rov(j)]/lamem(s)*2*pi;
            nv = [kbead kem kwater];
        end
        resx = []; resz = [];      
        for k = 1:length(rm)
            r = rm(k)/lamem(s)*2*pi;
            
            % coefficients of the dipole fields, as required by Spherical
            % x-orientation
            if r<rv(1)
                cx = [zeros(4*length(rv)-2,N); -nv(1)*sqrt(pi*(2*n+1)./n./(n+1)).*bj(n,nv(1)*r); ...
                        -i*nv(1)*sqrt(pi*(2*n+1)./n./(n+1)).*bjx(n,nv(1)*r)];
            else                
                cx = [-nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bh(n,nv(end)*r); ...
                        -i*nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bhx(n,nv(end)*r); zeros(4*length(rv)-2,N)];
            end 

            % z-orientation
            if r<rv(1)    
                cz = [zeros(4*length(rv)-1,N); i*sqrt(4*pi*(2*n+1)).*bj(n,nv(1)*r)/r];
            else
                cz = [zeros(1,N); i*sqrt(4*pi*(2*n+1)).*bh(n,nv(end)*r)/r; zeros(4*length(rv)-2,N)];        
            end

            if rm(k)<rov(j)-d | r>rv(end)
                coefx = Spherical(n,nv,rv,cx);
                coefz = Spherical(n,nv,rv,cz);
                
                % emission intensity normalized to polymer
                if r<rv(1)
                    resx(k,1) = 2*SphericalFlux(coefx(end-1,:),coefx(end,:),cx(end-1,:),cx(end,:),nv(1),rv(1));
                    resx(k,2) = 2*SphericalFlux(zeros(1,N),zeros(1,N),coefx(1,:),coefx(2,:),nv(end),rv(end));
                    resz(k,1) = SphericalFlux(coefz(end-1,:),coefz(end,:),cz(end-1,:),cz(end,:),nv(1),rv(1));
                    resz(k,2) = SphericalFlux(zeros(1,N),zeros(1,N),coefz(1,:),coefz(2,:),nv(end),rv(end));
                else
                    cxh = [-nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bj(n,nv(end)*r); ...
                        -i*nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bjx(n,nv(end)*r)];
                    resx(k,2) = 2*SphericalFlux(zeros(1,N),zeros(1,N),coefx(1,:)+cxh(1,:),coefx(2,:)+cxh(2,:),nv(end),inf);
                    resx(k,1) = resx(k,2) - 2*SphericalFlux(cx(1,:),cx(2,:),coefx(1,:),coefx(2,:),nv(end),rv(end));
                    czh = i*sqrt(4*pi*(2*n+1)).*bj(n,nv(end)*r)/r;        
                    resz(k,2) = SphericalFlux(zeros(1,N),zeros(1,N),zeros(1,N),czh+coefz(2,:),nv(end),inf);
                    resz(k,1) = resz(k,2) - SphericalFlux(cz(1,:),cz(2,:),coefz(1,:),coefz(2,:),nv(end),rv(end));                    
                end        
            else
                resx(k,1:2) = inf;
                resz(k,1:2) = inf;
            end
        end    
        qux(:,j,s) = resx(:,2)./resx(:,1); 
        emx(:,j,s) = resx(:,2)/(kwater/3);
        ltx(:,j,s) = kwater/3./resx(:,1);
        quz(:,j,s) = resz(:,2)./resz(:,1); 
        emz(:,j,s) = resz(:,2)/(kwater/3);
        ltz(:,j,s) = kwater/3./resz(:,1);
    end
end
save copperlifetime qux quz emx emz ltx ltz rov rm lamex lamem N kbead kwater d
 
c0 = 2*pi*i.^(n+1).*sqrt((2*n+1)/4/pi./n./(n+1));
nn = (n'.*(n+1)')*ones(1,length(rm));
for s=1:length(lamex)
    kex = copper(wavelength==lamex(s));
    nv = [kbead kex kwater];
    r = rm/lamex(s)*2*pi;
    for j = 1:length(rov)
        disp(['ex = ' num2str(lamex(s)) ' & j = ' num2str(j)]);
        if rov(j)<=d
            rv = rov(j)/lamex(s)*2*pi;
            c = [c0; c0; zeros(4*length(rv)-2,N)];
            coef = Spherical(n,nv(2:end),rv,c);
        else            
            rv = [rov(j)-d rov(j)]/lamex(s)*2*pi;
            c = [c0; c0; zeros(4*length(rv)-2,N)];
            coef = Spherical(n,nv,rv,c);
        end
        ind1 = rm<rov(j)-d; ind2 = r>rv(end); ind3 = ~(ind1 | ind2);        
        for k=1:N
            cm(k,ind1) = coef(end-1,k)*bj(k,nv(1)*r(ind1));
            cm(k,ind3) = 0;
            cm(k,ind2) = coef(1,k)*bh(k,nv(end)*r(ind2)) + c(1,k)*bj(k,nv(end)*r(ind2));
            cn(k,ind1) = coef(end,k)*bjx(k,nv(1)*r(ind1));
            cn(k,ind3) = 0;
            cn(k,ind2) = coef(2,k)*bhx(k,nv(end)*r(ind2)) + c(1,k)*bjx(k,nv(end)*r(ind2));
            cr(k,ind1) = coef(end,k)*bj(k,nv(1)*r(ind1))./r(ind1)/nv(1);
            cr(k,ind3) = 0;
            cr(k,ind2) = coef(2,k)*bh(k,nv(end)*r(ind2))./r(ind2)/nv(end) + c(1,k)*bj(k,nv(end)*r(ind2))./r(ind2)/nv(end);
        end
        intx(j,:,s) = sum(nn.*(abs(cm).^2 + abs(cn).^2))/4/pi;
        intz(j,:,s) = sum(abs(nn.*cr).^2)/2/pi;
    end
end
save copperlifetime intx intz -append
