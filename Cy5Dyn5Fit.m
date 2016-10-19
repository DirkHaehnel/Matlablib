function [err, c, zz] = Cy5Dyn5Fit(p, t, y);

% details see cy5dynamics.nb

kiso = p(3);
kbiso = p(4);
ksn = p(5);
kns = p(6);
ksniso = p(7);
knsiso = p(8);
mm = [-kiso-ksn, kns, kbiso, 0; ksn, -kns, 0, 0; kiso, 0, -kbiso-ksniso, knsiso; 0, 0, ksniso, -knsiso];
[ev,lambda] = eig(mm);
ev = ev';
lambda = diag(lambda)';
iv = inv(ev);

if length(t)<length(y)
    c = t;
    t = y;
    tmp = exp(t*lambda);
    shorttime = [((ones(size(t))*iv(1,:)).*tmp)*ev(:,1) ((ones(size(t))*iv(1,:)).*tmp)*ev(:,3) ...
            ((ones(size(t))*iv(3,:)).*tmp)*ev(:,1) ((ones(size(t))*iv(3,:)).*tmp)*ev(:,3)];
    col = ones(length(t),1);    
    z = ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,4)).*shorttime;
    for k=1:4
        err(:,k) = [col z(:,k)]*c(:,k);
    end
else
    t = t(:);
    tmp = exp(t*lambda);
    shorttime = [((ones(size(t))*iv(1,:)).*tmp)*ev(:,1) ((ones(size(t))*iv(1,:)).*tmp)*ev(:,3) ...
            ((ones(size(t))*iv(3,:)).*tmp)*ev(:,1) ((ones(size(t))*iv(3,:)).*tmp)*ev(:,3)];
    col = ones(length(t),1);
    z = ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,4)).*shorttime;
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
