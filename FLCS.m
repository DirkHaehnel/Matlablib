function [err, c, z] = FLCS(p, t, y, bld)

% Function [err, c, z] = FLCS(p, t, y) calculates the error err, amplitudes c, and fit curve z for the 
% FCS data y with the model offset + 1/(p(1)+t)/sqrt(p(2)+t)*(1 + exp(-t/p(3)) + exp(-t/p(4)) + exp(-t/p(5)))
% input y has to be a n-by-4 matrix

if length(t)<size(y,1)
    c = t;
    t = y(:);
    p = p(:)';
    z = [ones(length(t),1) ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,length(p)-1)).*[ones(length(t),1) exp(-t*(1./p(3:end)))]];
    err = z*c;
else
    t = t(:);
    p = p(:)';
    y = y(:,:);
    col = ones(size(t));
    z = [ones(length(t),1) ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,length(p)-2)).*[col exp(-t*(1./p(3:end-1)))]];
    zz = (p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t)).*exp(-t/p(end));
    if isreal(z)
        c(:,1) = lsqnonneg([z zz],y(:,1));
        c(:,2) = lsqnonneg([z -zz],y(:,2)); c(end,2) = -c(end,2);
        c(:,3) = lsqnonneg([z -zz],y(:,3)); c(end,3) = -c(end,3);
        c(:,4) = lsqnonneg([z zz],y(:,4));
        z = [z zz]*c;
        if nargin>3 && ~isempty(bld)
            semilogx(t,y./(col*mean(y)),'--',t,z./(col*mean(y))); drawnow;
        end
        err = sum(sum(y-z,2).^2./abs(sum(z,2)));
    else
        err = inf;
    end
end