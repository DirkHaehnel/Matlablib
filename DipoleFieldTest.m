% DipoleFieldTest

if 0
    lambda = 570;
    rhoa = 2*pi*(1:0.001:2)/lambda;
    zd = 0; za = 5*2*pi/lambda;
    n0 = 2;
    nd = 2;
    ni = 2;
    na = 2;
    n1 = 2;
    d0 = [];
    dd = 0;
    di = 0.02*pi/lambda;
    da = 0;
    d1 = [];

    [ex,ez,ex0,ez0] = DipoleField(zd,rhoa,za,n0,nd,ni,na,n1,d0,dd,di,da,d1);

    rr = sqrt((dd-zd + sum(di) + za)^2+rhoa.^2);
    ct = (dd-zd + sum(di) + za)./rr; st = rhoa./rr;
    exx(1,:) = (nd^2./rr + nd*1i./rr.^2 - 1./rr.^3 - 0.5*(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*st.^2).*exp(1i*rr)/nd^2;
    exx(2,:) = (-0.5*(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*st.^2).*exp(1i*rr)/nd^2;
    exx(3,:) = (-0.5*(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*st.*ct).*exp(1i*rr)/nd^2;
    ezz(1,:) = (nd^2./rr + nd*1i./rr.^2 - 1./rr.^3 -(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*ct.^2).*exp(1i*rr)/nd^2;
    ezz(2,:) = (-(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*st.*ct).*exp(1i*rr)/nd^2;
    plot(rhoa,abs(ezz).^2,'o',rhoa,abs(ez).^2)
    figure
    plot(rhoa,abs(exx).^2,'o',rhoa,abs(ex).^2)
end

if 1 % FRET over glass surface in water (neglecting donor lifetime)
    close all
    lambda = 570;
    rad = 2*pi/lambda*0.001*exp(log(3/0.001)*(0:0.01:1));
    phi = (0:90)/180*pi;
    zd = 0;
    n0 = 1.51;
    nd = 1.33;
    ni = 1.33;
    na = 1.33;
    n1 = 1.33;
    d0 = [];
    dd = zd;
    di = 0.2*pi/lambda;
    da = 0;
    d1 = [];
    
    ex = zeros(3,length(rad),length(phi));
    ex0 = ex;
    ez = zeros(2,length(rad),length(phi));
    ez0 = ez;
    
    for jr = 1:length(rad)
        for jp = 1:length(phi)
            [ex(:,jr,jp),ez(:,jr,jp),ex0(:,jr,jp),ez0(:,jr,jp)] = DipoleField(zd,rad(jr)*sin(phi(jp)),rad(jr)*cos(phi(jp)),n0,nd,ni,na,n1,d0,dd,di,da,d1);
        end
        jr
    end
    
    extotal = squeeze(abs(ex(1,:,:)).^2+0.5*abs(ex(2,:,:)).^2+0.5*abs(ex(3,:,:)).^2);
    ex0total = squeeze(abs(ex0(1,:,:)).^2+0.5*abs(ex0(2,:,:)).^2+0.5*abs(ex0(3,:,:)).^2);
    eztotal = squeeze(abs(ez(1,:,:)).^2+0.5*abs(ez(2,:,:)).^2);
    ez0total = squeeze(abs(ez0(1,:,:)).^2+0.5*abs(ez0(2,:,:)).^2);
    extotal = extotal + ex0total;
    eztotal = eztotal + ez0total;
    
    pcolor(rad'*sin(phi)/2/pi*lambda,rad'*cos(phi)/2/pi*lambda,extotal./ex0total); axis image; shading flat; colorbar
    figure
    pcolor(rad'*sin(phi)/2/pi*lambda,rad'*cos(phi)/2/pi*lambda,eztotal./ez0total); axis image; shading flat; colorbar
    figure
    pcolor(rad'*sin(phi)/2/pi*lambda,rad'*cos(phi)/2/pi*lambda,(2*extotal+eztotal)./(2*ex0total+ez0total)); axis image; shading flat; colorbar
    
    cameratoolbar; surf(rad'*sin(phi)/2/pi*lambda,rad'*cos(phi)/2/pi*lambda,(2*extotal+eztotal)./(2*ex0total+ez0total)); shading interp; view(115,16); alpha(0.8)
    cameratoolbar; surf(rad'*sin(phi)/2/pi*lambda,rad'*cos(phi)/2/pi*lambda,extotal./ex0total); camlight; shading interp; view(115,16); alpha(0.8)
    cameratoolbar; surf(rad'*sin(phi)/2/pi*lambda,rad'*cos(phi)/2/pi*lambda,eztotal./ez0total); camlight; shading interp; view(115,16); alpha(0.8)
    %plot(rhoa/2/pi*lambda,ez0total,rhoa/2/pi*lambda,eztotal,'o-')
    %plot(rhoa/2/pi*lambda,ex0total,rhoa/2/pi*lambda,extotal,'o-')
end

if 0 % FRET over gold and glass surface in water
    close all
    clear all
    lambda = 570;
    load metals
    nau = gold(wavelength==lambda);
    rhoa = 0;
    zd = 10/lambda*2*pi;
    zav = 2*pi*(0.2:0.01:10)/lambda;
    n0 = [1.51 nau];
    nd = 1.33;
    ni = 1.33;
    na = 1.33;
    n1 = 1.33;
    d0 = 25/lambda*2*pi;
    dd = zd;
    di = 0.2*pi/lambda;
    da = 0;
    d1 = [];

    for k=1:length(zav)
        [ex(:,k),ez(:,k),ex0(:,k),ez0(:,k)] = DipoleField(zd,rhoa,zav(k),n0,nd,ni,na,n1,d0,dd,di,zav(k),d1);
    end
    
    if 0
        semilogy(rhoa/2/pi*lambda,abs(ez0).^2,rhoa/2/pi*lambda,abs(ez+ez0).^2,'o-')
        figure
        semilogy(rhoa/2/pi*lambda,abs(ex0).^2,rhoa/2/pi*lambda,abs(ex+ex0).^2,'o-')
    else
        plot(zav/2/pi*lambda,abs(ex+ex0).^2./abs(ex0).^2,zav/2/pi*lambda,ones(size(rhoa)),':')
        figure
        plot(zav/2/pi*lambda,abs(ez+ez0).^2./abs(ez0).^2,zav/2/pi*lambda,ones(size(rhoa)),':')
    end
    %plot(rhoa/2/pi*lambda,sum(abs(ez0).^2),rhoa/2/pi*lambda,sum(abs(ez).^2),'o-')
    %plot(rhoa/2/pi*lambda,sum(abs(ex0).^2),rhoa/2/pi*lambda,sum(abs(ex).^2),'o-')
end

