function [err, c, z] = RiglerTriplet(p, t, y)

if length(t)<length(y)
    c = t;
    if size(y,1)==1
        t = y(:);
    else
        t = y;
    end
    p = p(:)';
    z = [ones(length(t),1) ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,length(p)-1)).*[ones(length(t),1) exp(-t*p(3:end))]];
    err = z*c;
else
    t = t(:);
    if size(y,1)==1
        y = y(:);
    end
    p = p(:)';
    len = (length(p)-2)/size(y,2);
    z0 = [ones(length(t),1) p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t)];
    for j=1:size(y,2)
        tmp = [z0 (z0(:,2)*ones(1,len)).*exp(-t*p(3+(j-1)*len:2+j*len))];
        c = lsqnonneg(tmp,y(:,j));
        z(:,j) = tmp*c;
    end
    semilogx(t,y,'--',t,z); drawnow;
    err = sum(sum((y-z).^2./abs(z)))
end
