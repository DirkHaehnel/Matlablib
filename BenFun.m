function [err, c, zz, z] = BenFun(p, t, y, pic)

if length(t)<length(y)
    c = t;
    t = y(:);
    p = p(:)';
    t = t(:)-min(t); p = 1./p(:)';
    col = ones(length(t),1);
    zz = diff((col*(1./p)).*exp(-t*p),1,2);
    zz = zz.*sign(zz);
	err = [col zz]*c;
else
    t = t(:)-min(t); p = 1./p(:)';
    y = y(:);
    col = ones(length(t),1);
    zz = diff((col*(1./p)).*exp(-t*p),1,2);
    zz = zz.*sign(zz);
    c = [col zz]\y;
	z = [col zz]*c;

    if nargin>3 && ~isempty(pic)
        %semilogy(t, y, 'o', t, z); drawnow
        plot(t, z, t, y, 'o'); drawnow
    end

    err = sum(sum((y-z).^2./abs(z)));
    %err = sum((y-z).^2);
    %ind = y>0;
    %err = sum(y(ind).*log(y(ind)./z(ind))-y(ind)+z(ind));
end

