function [err, c, z] = RiglerQD(p, t, y);

% Function [err, c, z] = RiglerQD(p, t, y) calculates the error err, amplitudes c, and fit curve z for the 
% FCS data y with the model offset + 1/(p(1)+t)/sqrt(p(2)+t)*(1 + t.^(-p(3)) + ...)

if length(t)<length(y)
    c = t;
    t = y(:);
    p = p(:)';
    z = [ones(length(t),1) ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,length(p)-1)).*[ones(length(t),1) exp(-log(t)*p(3:end))]];
    err = z*c;
else
    t = t(:);
    y = y(:);
    p = p(:)';
    z = [ones(length(t),1) ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,length(p)-1)).*[ones(length(t),1) exp(-log(t)*p(3:end))]];
    if isreal(z)
        c = lsqnonneg(z,y);
        z = z*c;
        semilogx(t,y,'--',t,z); drawnow;
        % err = sum((y-z).^2./abs(z))
        err = sum((y-z).^2)
    else
        err = inf;
    end
end
