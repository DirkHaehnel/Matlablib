function [err, c, z] = rigler(p, t, y, m);

% Function [err, c, z] = rigler(p, t, y) calculates the error err, amplitudes c, and fit curve z for the 
% FCS data y with the model offset + 1/(p(1)+t)/sqrt(p(2)+t)*(1 + exp(-t/p(3)) + ...)

if length(t)<length(y)
    c = t;
    t = y(:);
    p = p(:)';
    if nargin>3 & ~isempty(m)
        z = [ones(length(t),1) (p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))];
        for j=1:m-1
            z = [z (p(2+j)^(3/2)*p(1)*sqrt(p(2)))./(p(2+j)*p(1)+t)./sqrt(p(2+j)*p(2)+t)];
        end
        z = [z exp(-t*p(m+2:end))];
    else
        z = [ones(length(t),1) ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,length(p)-1)).*[ones(length(t),1) exp(-t*p(3:end))]];
    end
    err = z*c;
else
    t = t(:);
    y = y(:);
    p = p(:)';
    if nargin>3 & ~isempty(m)
        z = [ones(length(t),1) (p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))];
        for j=1:m-1
            z = [z (p(2+j)^(3/2)*p(1)*sqrt(p(2)))./(p(2+j)*p(1)+t)./sqrt(p(2+j)*p(2)+t)];
        end
        z = [z exp(-t*p(m+2:end))];
    else
        z = [ones(length(t),1) ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,length(p)-1)).*[ones(length(t),1) exp(-t*p(3:end))]];
    end
    if isreal(z)
        c = lsqnonneg(z,y);
        z = z*c;
        semilogx(t,y,'--',t,z); drawnow;
        err = sum((y-z).^2./abs(z))
        % err = sum((y-z).^2)
    else
        err = inf;
    end
end
