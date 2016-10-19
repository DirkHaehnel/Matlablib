function [err, c, z] = RiglerFlow(p, t, y, m, bld)

% Function [err, c, z] = RiglerFlow(p, t, y) calculates the error err, amplitudes c, and fit curve z for the 
% FCS data y with the model offset + 1./(p(1)+t)./sqrt(p(2)*p(1)+t).*exp(-p(3)*t.^2/(p(1)+t)).*(1 + exp(-t/p(4)) + ...)

if length(t)<length(y)
    c = t;
    t = y(:);
    p = p(:)';
    if nargin>3 && ~isempty(m)
        z = [ones(length(t),1) ((p(1)*sqrt(p(2)*p(1))./(p(1)+t)./sqrt(p(2)*p(1)+t).*exp(-p(3)*t.^2./(p(1)+t)))*ones(1,length(p)-m-1)).*[ones(length(t),1) exp(-t*p(3+m:end))]];
        for j=1:m-1
            z = [z (((p(2+j)^(3/2)*p(1)*sqrt(p(2)*p(1)))./(p(2+j)*p(1)+t)./sqrt(p(2+j)*p(2)*p(1)+t).*exp(-p(3)*t.^2./(p(1)+t)))*ones(1,length(p)-m-1)).*[ones(length(t),1) exp(-t*p(3+m:end))]];
        end
    else
        z = [ones(length(t),1) ((p(1)*sqrt(p(2)*p(1))./(p(1)+t)./sqrt(p(2)*p(1)+t).*exp(-p(3)*t.^2./(p(1)+t)))*ones(1,length(p)-2)).*[ones(length(t),1) exp(-t*p(4:end))]];
    end
    err = z*c;
else
    t = t(:);
    y = y(:);
    p = p(:)';
    if nargin>3 && ~isempty(m)
        z = [ones(length(t),1) ((p(1)*sqrt(p(2)*p(1))./(p(1)+t)./sqrt(p(2)*p(1)+t).*exp(-p(3)*t.^2./(p(1)+t)))*ones(1,length(p)-m-1)).*[ones(length(t),1) exp(-t*p(3+m:end))]];
        for j=1:m-1
            z = [z (((p(2+j)^(3/2)*p(1)*sqrt(p(2)*p(1)))./(p(2+j)*p(1)+t)./sqrt(p(2+j)*p(2)*p(1)+t).*exp(-p(3)*t.^2./(p(1)+t)))*ones(1,length(p)-m-1)).*[ones(length(t),1) exp(-t*p(3+m:end))]];
        end
    else
        z = [ones(length(t),1) ((p(1)*sqrt(p(2)*p(1))./(p(1)+t)./sqrt(p(2)*p(1)+t).*exp(-p(3)*t.^2./(p(1)+t)))*ones(1,length(p)-2)).*[ones(length(t),1) exp(-t*p(4:end))]];
    end
    if isreal(z)
        c = lsqnonneg(z,y);
        % c = z\y;
        z = z*c;
        if nargin>4 && ~isempty(bld)
            semilogx(t,y,'o',t,z, 'markersize', 2.5); drawnow;
        end
        %err = sum((y-z).^2./abs(z));
        err = sum((y-z).^2);
    else
        err = inf;
    end
end
