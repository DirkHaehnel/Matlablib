function [err, c, zz] = Cy5Dyn3Fit(p, t, y, p0);

% details see cy5dynamics.nb

kiso = p(1);
kbiso = p0(3)-p(1);
ksn = p(2);
kns = p0(4)-p(2);

dd = sqrt((kiso+kbiso+ksn+kns)^2-4*(kns*kiso+kbiso*(ksn+kns)));
lambda = -[kiso+kbiso+ksn+kns+dd, kiso+kbiso+ksn+kns-dd]/2;
ev(1,:) = [kbiso/kiso, kbiso/kiso*ksn/kns, 1];
ev(2,:) = [-(-kbiso + kiso + kns + ksn + dd)/(2*kiso), (-kbiso - kiso + kns + ksn + dd)/(2*kiso), 1];
ev(3,:) = [-(-kbiso + kiso + kns + ksn - dd)/(2*kiso), (-kbiso - kiso + kns + ksn - dd)/(2*kiso), 1];

iv = inv(ev);

if length(t)<length(y)
    c = t;
    t = y;
    tmp = [ones(size(t)) exp(t*lambda)];
    shorttime = [((ones(size(t))*iv(1,:)).*tmp)*ev(:,1) ((ones(size(t))*iv(1,:)).*tmp)*ev(:,3) ...
            ((ones(size(t))*iv(3,:)).*tmp)*ev(:,1) ((ones(size(t))*iv(3,:)).*tmp)*ev(:,3)];
    col = ones(length(t),1);    
    z = ((p0(1)*sqrt(p0(2))./(p0(1)+t)./sqrt(p0(2)+t))*ones(1,4)).*shorttime;
    for k=1:4
        err(:,k) = [col z(:,k)]*c(:,k);
    end
else
    t = t(:);
    tmp = [ones(size(t)) exp(t*lambda)];
    shorttime = [((ones(size(t))*iv(1,:)).*tmp)*ev(:,1) ((ones(size(t))*iv(1,:)).*tmp)*ev(:,3) ...
            ((ones(size(t))*iv(3,:)).*tmp)*ev(:,1) ((ones(size(t))*iv(3,:)).*tmp)*ev(:,3)];
    col = ones(length(t),1);
    z = ((p0(1)*sqrt(p0(2))./(p0(1)+t)./sqrt(p0(2)+t))*ones(1,4)).*shorttime;
    if isreal(z)
        for k=1:4
            c(:,k) = lsqnonneg([col z(:,k)],y(:,k));
            zz(:,k) = [col z(:,k)]*c(:,k);
        end
        semilogx(t,y,'--',t,zz); drawnow;
        err = sum(sum((y-zz).^2.))
    else
        err = inf;
    end
end
