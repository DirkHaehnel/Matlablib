function [z, zc] = BindingFCSXGlobal(p, delta, lam, zd, t, y1, y2, yc)

% p(1:2) - laser beam waist 
% p(3:4) - laser beam waist position
% p(5) - kp
% p(6) - km
% p(7) - diffusion coefficient

t = t(:); 
p = p(:);

if length(p)==5
    p = p([1:2 2 2 3:end]);
    p(3:4) = 0;
end

w1 = LaserBeamFit(p([1 3]),lam,zd).^2;
w2 = LaserBeamFit(p([2 4]),lam,zd).^2;
for k=1:length(zd)
    if p(5)<=0 | p(6)<0
        for k=1:length(zd)
            z1(:,1,k) = 1./(4*p(7)*t+w1(k));
            z2(:,1,k) = 1./(4*p(7)*t+w2(k));
            zc(:,1,k) = 1./(4*p(7)*t+w1(k)/2+w2(k)/2).*exp(-delta^2./(4*p(7)*t+w1(k)/2+w2(k)/2));
        end
    else
        mm = 100;
        q = (0.5:mm)/mm*5/max(p(1:2));
        kp = p(5)+p(6);
        km = p(5)-p(6);
        kk = km/kp^2*p(7);
        tmp = sqrt((p(7)*q.^2+km).^2+4*p(5)*p(6))/kp;
        om1 = 0.5*(p(7)*q.^2+kp*(1-tmp));
        om2 = 0.5*(p(7)*q.^2+kp*(1+tmp));
        for j=1:length(t)
            v = (q.*((1+tmp+kk*q.^2).*exp(-om1*t(j))-(1-tmp+kk*q.^2).*exp(-om2*t(j)))./(2*tmp));
            z1(j,1,k) = v*(exp(-w1(k)*q.^2/4))';
            z2(j,1,k) = v*(exp(-w2(k)*q.^2/4))';            
            zc(j,1,k) = v*(besselj(0,delta*q).*exp(-(w1(k)+w2(k))*q.^2/8))';
        end
    end
end

if nargin>5
    if isreal(z1) & isreal(z2)
        col = ones(length(t),1);
        for k=1:length(zd)
            c1 = lsqnonneg([col z1(:,:,k)],y1(:,k));
            ym1(:,k) = [col z1(:,:,k)]*c1;
            c2 = lsqnonneg([col z2(:,:,k)],y2(:,k));
            ym2(:,k) = [col z2(:,:,k)]*c2;
            ymc(:,k) = [col zc(:,:,k)]*sqrt(c1.*c2);
        end
        if length(zd)==1
            semilogx(t,[y1 y2 yc],'o',t,[ym1 ym2 ymc]); drawnow;
        else
            semilogx(t,y1./[col*ym1(end,:)],'o',t,ym1./[col*ym1(end,:)]); drawnow;
        end
        z = sum(sum(([y1 y2 yc]-[ym1 ym2 ymc]).^2./abs([ym1 ym2 ymc])));
        zc = [ym1 ym2 ymc];
    else 
        z=inf; zc = [];
    end
end
