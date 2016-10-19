function [err, c, z] = Rigler2f(p, t, y, yx);

% Function [err, c, z] = Rigler2D(p, t, y) calculates the error err, amplitudes c, and fit curve z for the 
% FCS data y with the model offset + 1/(p(1)+t)*(1 + exp(-t*p(2:end)) + ...)

if length(t)<size(y,1)
    c = t;
    t = y;
    p = p(:)';
    col = ones(length(t),1);    
    z1 = ((p(1)./(p(1)+t))*ones(1,length(p)-2)).*[ones(length(t),1) exp(-t*(1./p(4:end)))];
    z2 = ((p(2)./(p(2)+t))*ones(1,length(p)-2)).*[ones(length(t),1) exp(-t*(1./p(4:end)))];
    zc = ((0.5*(p(1)+p(2))./((p(1)+p(2))/2+t).*exp(-p(3)./((p(1)+p(2))/2+t)))*ones(1,length(p)-2)).*[ones(length(t),1) exp(-t*(1./p(4:end)))];
	err = [[col z1]*c(:,1) [col z2]*c(:,2) [col zc]*sqrt(c(:,1).*c(:,2))];
else
    t = t(:);
    p = p(:)';
    col = ones(length(t),1);
    z1 = ((p(1)./(p(1)+t))*ones(1,length(p)-2)).*[ones(length(t),1) exp(-t*(1./p(4:end)))];
    z2 = ((p(2)./(p(2)+t))*ones(1,length(p)-2)).*[ones(length(t),1) exp(-t*(1./p(4:end)))];
    zc = ((0.5*(p(1)+p(2))./((p(1)+p(2))/2+t).*exp(-p(3)./((p(1)+p(2))/2+t)))*ones(1,length(p)-2)).*[ones(length(t),1) exp(-t*(1./p(4:end)))];
    c(:,1) = lsqnonneg([col z1],y(:,1));
    c(:,2) = lsqnonneg([col z2],y(:,2));    
    z(:,1) = [col z1]*c(:,1);
    z(:,2) = [col z1]*c(:,2);
    z(:,3) = [col zc]*sqrt(c(:,1).*c(:,2));
    y(:,3) = yx;
    semilogx(t,y,'--',t,z); drawnow;
    err = sum(sum((y-z).^2./abs(z)));
    err = sum(sum((y-z).^2))
end

