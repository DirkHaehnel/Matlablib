function [err, c1, c2, cc] = MultiFocusCrossFit(p, t, y1, yc);

if length(t)<length(yc)
    c1 = t;
    c2 = y1;
    cc = y2;
    t = yc(:);
    p = p(:)';
    z = ones(length(t),1); zz = z;
    for j=1:(length(p)-2)/2
        z = [z ((p(2*j+1)*sqrt(p(2*j+2))./(p(2*j+1)+t)./sqrt(p(2*j+2)+t))*ones(1,2)).*[ones(length(t),1) exp(-t*p(2))]];
        zz = [zz ((p(2*j+1)*sqrt(p(2*j+2))./(p(2*j+1)+t)./sqrt(p(2*j+2)+t).*exp(-p(1)./(p(2*j+1)+t)))*ones(1,2)).*...
            [ones(length(t),1) exp(-t*p(2))]];
    end
    err = [z*c1 z*c2 zz*cc];  
else
    t = t(:);
    y1 = y1(:);
    yc = yc(:);
    p = p(:)';

%     z = 1./(p(2)+t)./sqrt(p(3)+t);
%     for j=1:3
%         z = [z exp(-0.4^2/4
	zz = y1.*exp(-p(1)./(p(2)+t));

    if length(p)>2
%         z = [z z.*exp(-t*p(3))];
        zz = [zz zz.*exp(-t*p(3))];
    end
%     z = [ones(size(z,1),1) z];
    zz = [ones(size(zz,1),1) zz];
    
    if isreal(zz)
%         c1 = lsqnonneg(z,y1);
%         c2 = lsqnonneg(z,y2);        
        cc = lsqnonneg(zz,yc);
%         z1 = z*c1;
%         z2 = z*c2;        
        zc = zz*cc;  
%         semilogx(t,y1,'--',t,z1,t,y2,'--',t,z2,t,yc,'--',t,zc); 
%         drawnow;
%         err = sum((y1-z1).^2.)+sum((y2-z2).^2.)+sum((yc-zc).^2.)
        semilogx(t,y1,'--',t,yc,'--',t,zc); 
        drawnow;
        err = sum((yc-zc).^2.)
    else
        err = inf;
    end
end
