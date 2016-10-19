close all

lamem = 0.67;
ns = 1.33;
n1 = 4;
n2 = 1.5;

theta = pi/2000:pi/1000:pi/2;
wave = 2*pi/lamem;
mm = 20;

if 0 % finding the bandgap
    w1 = sqrt(n1^2-ns^2*sin(theta).^2);
    w2 = sqrt(n2^2-ns^2*sin(theta).^2);
    [d1, d2] = meshgrid(wave*(0.01:0.001:0.2),wave*(0.01:0.001:0.2));
    j = 1;
    tst1 = double(abs(cos(w1(j)*d1).*cos(w2(j)*d2)-0.5*(n2^2/n1^2*w1(j)/w2(j)+n1^2/n2^2*w2(j)/w1(j)).*sin(w1(j)*d1).*sin(w2(j)*d2))>=1);
    tst2 = double(abs(cos(w1(j)*d1).*cos(w2(j)*d2)-0.5*(w1(j)/w2(j)+w2(j)/w1(j)).*sin(w1(j)*d1).*sin(w2(j)*d2))>1);
    for j=2:length(w1)
        tst1 = tst1.*double(abs(cos(w1(j)*d1).*cos(w2(j)*d2)-0.5*(n2^2/n1^2*w1(j)/w2(j)+n1^2/n2^2*w2(j)/w1(j)).*sin(w1(j)*d1).*sin(w2(j)*d2))>=1);
        tst2 = tst2.*double(abs(cos(w1(j)*d1).*cos(w2(j)*d2)-0.5*(w1(j)/w2(j)+w2(j)/w1(j)).*sin(w1(j)*d1).*sin(w2(j)*d2))>1);
    end
    pcolor(d1/wave,d2/wave,tst1+tst2);
    axis image; shading flat
end

if 0 % explicit reflectivity calculation
    cnt = 1;
    for fac = 0.05:0.0001:0.08
        d1 = fac*wave;
        d2 = fac*wave;

        if 1
            nv = [n1 repmat([n1 n2],1,mm) ns];
            dv = repmat([d1 d2],1,mm);
            [rp, rs] = Fresnel(ns*cos(theta), fliplr(nv), fliplr(dv));
            plot(theta/pi*180,abs(rp).^2,theta/pi*180,abs(rs).^2);
            axis([0 90 0 1]);
        end
        if 0 % dipole emission
            z = 0;
            nv = [n1 repmat([n1 n2],1,mm)];
            dv = repmat([d1 d2],1,mm);
            [vd,pcd,psd] = DipoleL(theta,z,nv,ns,ns,dv,z,[]);
            [vu,pcu,psu] = DipoleL(theta,0,ns,ns,fliplr(nv),[],z,fliplr(dv));
            subplot(121);
            mpolar([theta fliplr(pi-theta)],[abs(vd).^2; flipud(abs(vu).^2)]')
            title(['\itd\rm_1 = ' mnum2str(d1/wave,1,3)])
            subplot(122);
            mpolar([theta fliplr(pi-theta)],0.5*[abs(pcd).^2+abs(psd).^2; flipud(abs(pcu).^2+abs(psu).^2)]')
            title(['\itd\rm_2 = ' mnum2str(d2/wave,1,3)])
        end
        %mov(cnt) = getframe; cnt = cnt + 1;
        drawnow
    end
end

if 1 % search for optimal emission properties
    d1v = 0.05:0.001:0.1;
    d2v = 0.05:0.001:0.1;
    z = 0;
    nv = [n1 repmat([n1 n2],1,mm)];
    ss = sin(theta)*diff(theta(1:2));
    emivu = zeros(length(d2v),length(d1v));
    emivd = emivu; emivt = emivu;
    emipu = emivu; emipd = emivu; emipt = emivu;
    for j=1:length(d2v)
        d2 = d2v(j)*wave;
        for k=1:length(d1v)
            d1 = d1v(k)*wave;
            dv = repmat([d1 d2],1,mm);
            [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu,qv,qp] = LifetimeL(z,nv,ns,ns,dv,z,[]);
            emivu(j,k) = qvu;
            emivd(j,k) = qvd;            
            emipu(j,k) = qpu;
            emipd(j,k) = qpd;
            emivt(j,k) = qvu+qvd+qv;
            emipt(j,k) = qpu+qpd+qp;            
            [j k]
        end
    end
end
