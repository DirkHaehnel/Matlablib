function [err, c, z] = FLCS1(p, t, y, a, bld);

% Function [err, c, z] = FLCS1(p, t, y) calculates the error err, amplitudes c, and fit curve z for the 
% FCS data y with the model offset + 1/(p(1)+t)/sqrt(p(2)+t)*(1 + exp(-t/p(3)) + exp(-t/p(4)) + exp(-t/p(5)))
% input y has to be a n-by-4 matrix

if length(t)<size(y,1)
    c = t;
    t = y(:);
    p = p(:)';
    z = [ones(length(t),1) ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,length(p)-1)).*[col exp(-t*(1./p(3:4)))]];
    err = z*c;
else
    t = t(:);
    p = p(:)';
    y = y(:,:);
    col = ones(size(t));
    for j=1:4
        zz = [ones(length(t),1) ((a(1)*sqrt(a(2))*p(j)^(3/2)./(a(1)*p(j)+t)./sqrt(a(2)*p(j)+t))*ones(1,length(p)-3)).*[col exp(-t*(1./p(5:end)))]];
        c(:,j) = zz\y(:,j);
        z(:,j) = zz*c(:,j);
    end
    if nargin>4 & ~isempty(bld)
        semilogx(t,y./(col*mean(y)),'--',t,z./(col*mean(y))); drawnow;
    end
    err = sum(sum((y-z).^2./abs(z)));
end