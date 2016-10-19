function [z, zc] = BindingFCSGlobal(p, lam, zd, t, y1)

% p(1) - laser beam waist 
% p(2) - laser beam waist position
% p(3) - kp
% p(4) - km
% p(5) - diffusion coefficient

if length(p)==3
    p = p([1 1 1:end]);
    p(1) = 1;
    p(2) = 0;
end

t = t(:); 
p = p(:);

ww = LaserBeamFit(p(1:2),lam,zd).^2;
for k=1:length(zd)
    if p(3)<=0 | p(4)<0
        for k=1:length(zd)
            z(:,1,k) = 1./(4*p(5)*t+ww(k));
        end
    else
        mm = 100;
        q = (0.5:mm)/mm*5/p(1);
        kp = p(3)+p(4);
        km = p(3)-p(4);
        kk = km/kp^2*p(5);
        tmp = sqrt((p(5)*q.^2+km).^2+4*p(3)*p(4))/kp;
        om1 = 0.5*(p(5)*q.^2+kp*(1-tmp));
        om2 = 0.5*(p(5)*q.^2+kp*(1+tmp));
        for j=1:length(t)
            z(j,1,k) = (q.*((1+tmp+kk*q.^2).*exp(-om1*t(j))-(1-tmp+kk*q.^2).*exp(-om2*t(j)))./(2*tmp))*(exp(-ww(k)*q.^2/4))';
        end
    end
end
if length(p)==6
    for k=1:length(zd)
        z(:,2,k) = 1./(4*p(6)*t+ww(k));
    end
end

if nargin>4
    if isreal(z)
        col = ones(length(t),1);
        for k=1:length(zd)
            c1 = lsqnonneg([col z(:,:,k)],y1(:,k));
            ym1(:,k) = [col z(:,:,k)]*c1;
        end
        semilogx(t,y1./[col*ym1(end,:)],'o',t,ym1./[col*ym1(end,:)]); drawnow;
        z = sum(sum((y1-ym1).^2./abs(ym1)));
        zc = ym1;
    else 
        z=inf; zc = [];
    end
end
