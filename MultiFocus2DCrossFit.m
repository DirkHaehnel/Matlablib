function [err, c1, c2, cc] = MultiFocusCrossFit(p, t, y1, y2, yc);

if length(t)<length(yc)
    c1 = t;
    c2 = y1;
    cc = y2;
    t = yc(:);
    p = p(:)';
    z1 = [ones(length(t),1) ((p(1)./(p(1)+t))*ones(1,2)).*[ones(length(t),1) exp(-t*p(2))]];
    z2 = [ones(length(t),1) ((p(3)./(p(3)+t))*ones(1,2)).*[ones(length(t),1) exp(-t*p(4))]];
    zc = [ones(length(t),1) (((p(1)+p(3))./(p(1)+p(3)+2*t).*exp(-2*p(5)./(p(1)+p(3)+2*t)))*ones(1,2)).*...
            [ones(length(t),1) exp(-t*p(6))]];
    err = [z1*c1 z2*c2 zc*cc];  
else
    t = t(:);
    y1 = y1(:);
    y2 = y2(:);
    yc = yc(:);
    p = p(:)';
    z1 = [ones(length(t),1) ((p(1)./(p(1)+t))*ones(1,2)).*[ones(length(t),1) exp(-t*p(2))]];
    z2 = [ones(length(t),1) ((p(3)./(p(3)+t))*ones(1,2)).*[ones(length(t),1) exp(-t*p(4))]];
    zc = [ones(length(t),1) (((p(1)+p(3))./(p(1)+p(3)+2*t).*exp(-2*p(5)./(p(1)+p(3)+2*t)))*ones(1,2)).*...
            [ones(length(t),1) exp(-t*p(6))]];
    if isreal(z1) & isreal(z2)
        c1 = lsqnonneg(z1,y1);
        c2 = lsqnonneg(z2,y2);        
        cc = lsqnonneg(zc,yc);
        z1 = z1*c1;
        z2 = z2*c2;        
        zc = zc*cc;  
        semilogx(t,y1,'--',t,z1,t,y2,'--',t,z2,t,yc,'--',t,zc); 
        drawnow;
        err = sum((y1-z1).^2.)+sum((y2-z2).^2.)+sum((yc-zc).^2.);
    else
        err = inf;
    end
end
