function [err, c, z] = ResactFun(p, t, y, temp, bld);

% Function [err, c, z] = rigler(p, t, y) calculates the error err, amplitudes c, and fit curve z for the 
% FCS data y with the model offset + 1/(p(1)+t)/sqrt(p(2)+t)*(1 + exp(-t/p(3)) + ...)

if length(t)<length(y)
    c = t;
    t = y;
    p = p(:)';
    if nargin>3 & ~isempty(m)
        z = [ones(length(t),1) ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,length(p)-m)).*[ones(length(t),1) exp(-t*p(2+m:end))]];
        for j=1:m-1
            z = [z (((p(2+j)^(3/2)*p(1)*sqrt(p(2)))./(p(2+j)*p(1)+t)./sqrt(p(2+j)*p(2)+t))*ones(1,length(p)-m)).*[ones(length(t),1) exp(-t*p(2+m:end))]];
        end
    else
        z = [ones(length(t),1) ((p(1)*sqrt(p(2))./(p(1)+t)./sqrt(p(2)+t))*ones(1,length(p)-1)).*[ones(length(t),1) exp(-t*p(3:end))]];
    end
    err = z*c;
else
    t = t(:);
    p = p(:)';
    col = ones(length(t),1);
    kp = p(end-3)*exp(-p(end-1)./temp);
    km = p(end-2)*exp(-p(end)./temp);
    p = p(1:end-4);
    for j=1:size(y,2)
        z(:,j) = p(2*j-1)*sqrt(p(2*j))./(p(2*j-1)+t)./sqrt(p(2*j)+t).*(km(j)+kp(j)*exp(-t*(km(j)+kp(j))))/(km(j)+kp(j));
        c(:,j) = lsqnonneg([col z(:,j)],y(:,j));
        z(:,j) = [col z(:,j)]*c(:,j);
    end
    if nargin>4 & ~isempty(bld)
        semilogx(t,y,'o',t,z, 'markersize', 2.5); drawnow;
    end
    err = sum(sum((y-z).^2./abs(z)));
end

