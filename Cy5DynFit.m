function [err, c, zz] = Cy5DynFit(p, t, y);

% details see cy5dynamics.nb

if length(p)==8
    ksiso = 1/p(3);
    kniso = 1/p(4);
    ksph = 1/p(5);
    knph = 1/p(6);
    ksn = 1/p(7);
    kns = 1/p(8);
    M = [-ksiso-ksn, kns, ksph, 0;...
            ksn, -kniso-kns, 0, knph;...
            ksiso, 0, -ksn-ksph, kns;...
            0, kniso, ksn, -knph-kns];
else
    ksiso = 1/p(3);
    ksph = 1/p(4);
    ksn = 1/p(5);
    kns = 1/p(6);
    M = [-ksiso-ksn, kns, ksph, 0;...
            ksn, -ksiso-kns, 0, ksph;...
            ksiso, 0, -ksn-ksph, kns;...
            0, ksiso, ksn, -ksph-kns];
end

[ev,lambda] = eig(M);
lambda = diag(lambda)';
iv = inv(ev)';
ev = ev';

if length(t)<length(y)
    c = t;
    t = y;
    col = ones(length(t),1);    
    tmp = exp(t*lambda);
    shorttime = [((ones(size(t))*iv(1,:)).*tmp)*ev(:,1) ((ones(size(t))*iv(2,:)).*tmp)*ev(:,1) ...
            ((ones(size(t))*iv(1,:)).*tmp)*ev(:,2) ((ones(size(t))*iv(2,:)).*tmp)*ev(:,2)];
    z = ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,4)).*shorttime;
    for k=1:4
        err(:,k) = [col z(:,k)]*c(:,k);
    end
else
    t = t(:);
    y = y(:,:); 
    col = ones(length(t),1);    
    tmp = exp(t*lambda);
    shorttime = [((col*iv(1,:)).*tmp)*ev(:,1) ((col*iv(2,:)).*tmp)*ev(:,1) ...
            ((col*iv(1,:)).*tmp)*ev(:,2) ((col*iv(2,:)).*tmp)*ev(:,2)];
    if length(p)>8
        shorttime = [...
                shorttime(:,1) + p(9)*shorttime(:,3) + p(9)*shorttime(:,2) + p(9)^2*shorttime(:,4), ...
                shorttime(:,2) + p(9)*shorttime(:,4) + p(10)*shorttime(:,1) + p(9)*p(10)*shorttime(:,3), ...
                shorttime(:,3) + p(9)*shorttime(:,4) + p(10)*shorttime(:,1) + p(9)*p(10)*shorttime(:,2), ...
                shorttime(:,4) + p(10)*shorttime(:,3) + p(10)*shorttime(:,2) + p(10)^2*shorttime(:,1)];
    end
    col = ones(length(t),1);
    z = ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,4)).*shorttime;
    if isreal(z)
        for k=1:4
            c(:,k) = lsqnonneg([col z(:,k)],y(:,k));
            zz(:,k) = [col z(:,k)]*c(:,k);
        end
        semilogx(t,y./(col*mean(y)),'--',t,zz./(col*mean(y))); drawnow;
        %semilogx(t,y(:,4),'--',t,zz(:,4)); drawnow;
        err = prod(sum((y-zz).^2.))
    else
        err = inf;
    end
end

