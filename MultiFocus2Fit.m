function [err, c, z] = MultiFocus2Fit(p, t, y);

if length(t)<length(y)
    c = t;
    t = y(:);
    p = p(:)';
    z = [ones(length(t),1) ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t).*exp(-p(4)./(p(3)+t)))*ones(1,2))];
    err = z*c;  
else
    t = t(:);
    y = y(:);
    p = p(:)';
    z = [ones(length(t),1) 1./(p(1)+t).*exp(-p(3)./(p(2)+t))];
    if isreal(z)
        c = lsqnonneg(z,y);
        z = z*c;
        semilogx(t,y,'--',t,z); 
        drawnow;
        err = sum((y-z).^2)
    else
        err = inf;
    end
end
