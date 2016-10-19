% program for the calculation of the change in extinction

%%
close all
clear all
load metals
warning off

bj = inline('besselj(n+0.5,x).*sqrt(pi./x/2)','n','x');
bh = inline('besselh(n+0.5,x).*sqrt(pi./x/2)','n','x');
bjx = inline('sqrt(pi/8./x).*(besselj(n+0.5,x)./x + besselj(n-0.5,x) - besselj(n+1.5,x))','n','x');
bhx = inline('sqrt(pi/8./x).*(besselh(n+0.5,x)./x + besselh(n-0.5,x) - besselh(n+1.5,x))','n','x');

%%

outer_radius_vec = 2.5:0.2:30;
kbead = 1.5;
kwater = 1.0;
d = 5;
lambda = 500;
N = 30; n = 1:N;
rm = (0.5:200)/200*40;

kem = silver(wavelength==lambda);
kex = kem;

for outer_radius_cnt=1:length(outer_radius_vec)
    
    rov = outer_radius_vec(outer_radius_cnt);
    
    if rov<=d
        rv = rov/lambda*2*pi;
        nv = [kem kwater];
    else
        rv = [rov-d rov]/lambda*2*pi;
        nv = [kbead kem kwater];
    end
    resx = []; resz = [];
    for k = 1:length(rm)
        r = rm(k)/lambda*2*pi;
        
        % coefficients of the dipole fields, as required by Spherical
        % x-orientation
        if r<rv(1)
            cx = [zeros(4*length(rv)-2,N); -nv(1)*sqrt(pi*(2*n+1)./n./(n+1)).*bj(n,nv(1)*r); ...
                -1i*nv(1)*sqrt(pi*(2*n+1)./n./(n+1)).*bjx(n,nv(1)*r)];
        else
            cx = [-nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bh(n,nv(end)*r); ...
                -1i*nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bhx(n,nv(end)*r); zeros(4*length(rv)-2,N)];
        end
        
        % z-orientation
        if r<rv(1)
            cz = [zeros(4*length(rv)-1,N); 1i*sqrt(4*pi*(2*n+1)).*bj(n,nv(1)*r)/r];
        else
            cz = [zeros(1,N); 1i*sqrt(4*pi*(2*n+1)).*bh(n,nv(end)*r)/r; zeros(4*length(rv)-2,N)];
        end
        
        if rm(k)<rov-d || r>rv(end)
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
                    -1i*nv(end)*sqrt(pi*(2*n+1)./n./(n+1)).*bjx(n,nv(end)*r)];
                resx(k,2) = 2*SphericalFlux(zeros(1,N),zeros(1,N),coefx(1,:)+cxh(1,:),coefx(2,:)+cxh(2,:),nv(end),inf);
                resx(k,1) = resx(k,2) - 2*SphericalFlux(cx(1,:),cx(2,:),coefx(1,:),coefx(2,:),nv(end),rv(end));
                czh = 1i*sqrt(4*pi*(2*n+1)).*bj(n,nv(end)*r)/r;
                resz(k,2) = SphericalFlux(zeros(1,N),zeros(1,N),zeros(1,N),czh+coefz(2,:),nv(end),inf);
                resz(k,1) = resz(k,2) - SphericalFlux(cz(1,:),cz(2,:),coefz(1,:),coefz(2,:),nv(end),rv(end));
            end
        else
            resx(k,1:2) = inf;
            resz(k,1:2) = inf;
        end
    end
    emx(:,outer_radius_cnt) = resx(:,1);
    emz(:,outer_radius_cnt) = resz(:,1);
    infx(:,outer_radius_cnt) = resx(:,2);
    infz(:,outer_radius_cnt) = resz(:,2);

    
    %%
    
    c0 = 2*pi*1i.^(n+1).*sqrt((2*n+1)/4/pi./n./(n+1));
    nn = (n'.*(n+1)')*ones(1,length(rm));
    r = rm/lambda*2*pi;
    if rov<=d
        rv = rov/lambda*2*pi;
        c = [c0; c0; zeros(4*length(rv)-2,N)];
        coef = Spherical(n,nv(end-1:end),rv,c);
    else
        rv = [rov-d rov]/lambda*2*pi;
        c = [c0; c0; zeros(4*length(rv)-2,N)];
        coef = Spherical(n,nv,rv,c);
    end
    d0 = [sqrt(imag(conj(coef(1,:).*bh(n,nv(end)*rv(end)) + c(1,:).*bj(n,nv(end)*rv(end))).*...
        (coef(1,:).*bhx(n,nv(end)*rv(end)) + c(1,:).*bjx(n,nv(end)*rv(end)))))*nv(end)*rv(end); ...
        sqrt(imag(conj(coef(2,:).*bh(n,nv(end)*rv(end)) + c(2,:).*bj(n,nv(end)*rv(end))).*...
        (coef(2,:).*bhx(n,nv(end)*rv(end)) + c(2,:).*bjx(n,nv(end)*rv(end)))))*nv(end)*rv(end)];
    doef = zeros(8,N);
    if rov>d
        rr = rv(1)+(0.5:200)/200*(rv(2)-rv(1));
        for k=1:N
            fj = sum(bj(k,nv(2)*rr))*mean(diff(rr));
            fh = sum(bh(k,nv(2)*rr))*mean(diff(rr));
            fjx = sum(bjx(k,nv(2)*rr))*mean(diff(rr));
            fhx = sum(bhx(k,nv(2)*rr))*mean(diff(rr));
            M = [fj*bh(k,nv(2)*rv(2)), 0, bh(k,nv(2)*rv(2)), 0, bj(k,nv(2)*rv(2)), 0, 0, 0; ...
                0, fjx*bhx(k,nv(2)*rv(2)), 0, bhx(k,nv(2)*rv(2)), 0, bjx(k,nv(2)*rv(2)), 0, 0; ...
                fj*nv(2)*bhx(k,nv(2)*rv(2)), 0, nv(2)*bhx(k,nv(2)*rv(2)), 0, nv(2)*bjx(k,nv(2)*rv(2)), 0, 0, 0; ...
                0, fjx*nv(2)*bh(k,nv(2)*rv(2)), 0, nv(2)*bh(k,nv(2)*rv(2)), 0, nv(2)*bj(k,nv(2)*rv(2)), 0, 0;
                fh*bj(k,nv(2)*rv(1)), 0, bh(k,nv(2)*rv(1)), 0, bj(k,nv(2)*rv(1)), 0, -bj(k,nv(1)*rv(1)), 0; ...
                0, fhx*bjx(k,nv(2)*rv(1)), 0, bhx(k,nv(2)*rv(1)), 0, bjx(k,nv(2)*rv(1)), 0, -bjx(k,nv(1)*rv(1)); ...
                fh*nv(2)*bjx(k,nv(2)*rv(1)), 0, nv(2)*bhx(k,nv(2)*rv(1)), 0, nv(2)*bjx(k,nv(2)*rv(1)), 0, -nv(1)*bjx(k,nv(1)*rv(1)), 0; ...
                0, fhx*nv(2)*bj(k,nv(2)*rv(1)), 0, nv(2)*bh(k,nv(2)*rv(1)), 0, nv(2)*bj(k,nv(2)*rv(1)), 0, -nv(1)*bj(k,nv(1)*rv(1))];
            doef(:,k) = M\[d0(1,k)*bh(k,nv(3)*rv(2)); d0(2,k)*bhx(k,nv(3)*rv(2)); nv(3)*d0(1,k)*bhx(k,nv(3)*rv(2)); nv(3)*d0(2,k)*bh(k,nv(3)*rv(2)); zeros(4,1)];
        end
    end
    
    ind1 = rm<rov-d; ind2 = r>rv(end);
    cm = zeros(N,numel(r)); cn = cm; cr = cm;
    dm = cm; dn = cm; dr = cm;
    for k=1:N
        cm(k,ind1) = coef(end-1,k)*bj(k,nv(1)*r(ind1));
        cm(k,ind2) = coef(1,k)*bh(k,nv(end)*r(ind2)) + c(1,k)*bj(k,nv(end)*r(ind2));
        cn(k,ind1) = coef(end,k)*bjx(k,nv(1)*r(ind1));
        cn(k,ind2) = coef(2,k)*bhx(k,nv(end)*r(ind2)) + c(1,k)*bjx(k,nv(end)*r(ind2));
        cr(k,ind1) = coef(end,k)*bj(k,nv(1)*r(ind1))./r(ind1)/nv(1);
        cr(k,ind2) = coef(2,k)*bh(k,nv(end)*r(ind2))./r(ind2)/nv(end) + c(1,k)*bj(k,nv(end)*r(ind2))./r(ind2)/nv(end);
        
        dm(k,ind1) = doef(7,k)*bj(k,nv(1)*r(ind1));
        dn(k,ind1) = doef(8,k)*bjx(k,nv(1)*r(ind1));
        dr(k,ind1) = doef(8,k)*bj(k,nv(1)*r(ind1))./nv(1)./r(ind1);
        dm(k,ind2) = d0(1,k)*bh(k,nv(end)*r(ind2));
        dn(k,ind2) = d0(2,k)*bhx(k,nv(end)*r(ind2));
        dr(k,ind2) = d0(2,k)*bh(k,nv(end)*r(ind2))./nv(end)./r(ind2);
    end
    intx(:,outer_radius_cnt) = sum(nn.*(abs(cm).^2 + abs(cn).^2))'/4/pi*kwater;
    intz(:,outer_radius_cnt) = sum(abs(nn.*cr).^2)'/2/pi*kwater;
    flx(:,outer_radius_cnt) = sum(nn.*(abs(dm).^2 + abs(dn).^2))'/4/pi*kwater;
    flz(:,outer_radius_cnt) = sum(abs(nn.*dr).^2)'/2/pi*kwater;
    
    close all
    plot(rm,resx(:,1),rm,intx(:,outer_radius_cnt)+flx(:,outer_radius_cnt),'o')
    figure
    plot(rm,resz(:,1),rm,intz(:,outer_radius_cnt)+flz(:,outer_radius_cnt),'o')
    
    disp(['outer radius = ' mnum2str(rov,2,1)])
end

warning on

save CavityAbsorption emx emz infx infz intx intz flx flz rm outer_radius_vec d kbead kwater kem lambda N

