function [err, c, z] = RiglerAffine(p, t, y, bld);

% Function [err, c, z] = rigler(p, t, y) calculates the error err, amplitudes c, and fit curve z for the 
% FCS data y with the model offset + 1/(p(1)+t)/sqrt(p(2)+t)*(1 + exp(-t/p(3)) + ...)

m = size(y,2);
n = (length(p) - m - 1)/m;
if length(t)<length(y)
    c = t;
    t = y;
    p = p(:)';
    col = ones(length(t),1);
    zz = ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,n+1)).*[ones(length(t),1) exp(-t*p(2+m:m:end))];
    for j=1:m-1
        zz = [zz (((p(2+j)^(3/2)*p(1)*sqrt(p(2)))./(p(2+j)*p(1)+t)./sqrt(p(2+j)*p(2)+t))*ones(1,n+1)).*[ones(length(t),1) exp(-t*p(2+m+j:m:end))]];
    end
    for j=1:m
        % c = lsqnonneg(z,y);
        err(:,j) = [col zz(:,1+(j-1)*(n+1):j*(n+1))]*c(:,j);
    end
else
    t = t(:);
    p = p(:)';
    col = ones(length(t),1);
    zz = ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,n+1)).*[ones(length(t),1) exp(-t*p(2+m:m:end))];
    for j=1:m-1
        zz = [zz (((p(2+j)^(3/2)*p(1)*sqrt(p(2)))./(p(2+j)*p(1)+t)./sqrt(p(2+j)*p(2)+t))*ones(1,n+1)).*[ones(length(t),1) exp(-t*p(2+m+j:m:end))]];
    end
    if isreal(zz)
        for j=1:m
            % c = lsqnonneg(z,y);
            c(:,j) = [col zz(:,1+(j-1)*(n+1):j*(n+1))]\y(:,j);
            z(:,j) = [col zz(:,1+(j-1)*(n+1):j*(n+1))]*c(:,j);
        end
        if nargin>3 & ~isempty(bld)
            semilogx(t,y,'o',t,z, 'markersize', 2.5); drawnow;
        end
        % err = sum(sum((y-z).^2./abs(z)));
        err = sum(sum((y-z).^2));
    else
        err = inf;
    end
end
