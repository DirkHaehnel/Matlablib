function [err, c, z] = BleachScan(p, t, x);

t = t(:); 
x = x(:);
tt = (3*p(1)*(-24.5:25)/25);
[tmp,tt] = meshgrid(t,tt);
z = exp(-t.^2/2/p(1)^2).*sum(exp(-(tt+tmp).*tt/p(1)^2).*exp(-p(2)*(erf((tmp+tt)/sqrt(2)/p(1))-erf(tt/sqrt(2)/p(1)))))';
c = [ones(length(x),1) z]\x;
z = [ones(length(x),1) z]*c;

semilogx(t, z, t, x, 'o'); drawnow

err = sum((x-z).^2./abs(z))

