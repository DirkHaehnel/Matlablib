function [err, c, z] = Rigler2D(p, t, y, yx)

% Function [err, c, z] = Rigler2D(p, t, y) calculates the error err, amplitudes c, and fit curve z for the 
% FCS data y with the model offset + 1/(p(1)+t)*(1 + exp(-t*p(2:end)) + ...)

if length(t)<size(y,1)
    c = t;
    t = y;
    p = p(:)';
    if nargin>3
        z = [ones(length(t),1) ((p(1)./(p(1)+t))*ones(1,length(p)-1)).*[ones(length(t),1) exp(-t*(1./p(3:end)))]];
        zc = z; zc(:,2).*exp(-p(2)./(p(1)+t));        
        err = [z*c zc*sqrt(c(:,1).*c(:,2))]
    else
        z = [ones(length(t),1) ((p(1)./(p(1)+t))*ones(1,length(p))).*[ones(length(t),1) exp(-t*(1./p(2:end)))]];
        err = z*c;
    end
else
    t = t(:);
    p = p(:)';
    if nargin>3
        zz(:,:,1) = [ones(length(t),1) ((p(1)./(p(1)+t))*ones(1,length(p)-2)).*[ones(length(t),1) exp(-t*(1./p(4:end)))]];
        zz(:,:,2) = [ones(length(t),1) ((p(2)./(p(2)+t))*ones(1,length(p)-2)).*[ones(length(t),1) exp(-t*(1./p(4:end)))]];        
        zc = [ones(length(t),1) ((sqrt(p(1)*p(2))./sqrt((p(1)+t).*(p(2)+t)).*exp(-p(3)./(p(1)+t)))*ones(1,length(p)-2)).*[ones(length(t),1) exp(-t*(1./p(4:end)))]];
    else
        zz = [ones(length(t),1) ((p(1)./(p(1)+t))*ones(1,length(p))).*[ones(length(t),1) exp(-t*(1./p(2:end)))]];
    end
    if isreal(zz)
        for j=1:size(y,2)
            c(:,j) = lsqnonneg(zz(:,:,j),y(:,j));
            z(:,j) = zz(:,:,j)*c(:,j);
        end
        if nargin>3
            z(:,3) = zc*sqrt(c(:,1).*c(:,2));
            y(:,3) = yx;
        end
        semilogx(t,z,t,y,'o'); drawnow;
        err = sum(sum((y-z).^2./abs(z)))
        % err = sum((y-z).^2)
    else
        err = inf;
    end
end

