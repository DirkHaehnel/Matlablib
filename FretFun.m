function [err, c, zz, z] = FretFun(p, t, y, tau, pic)

if length(t)<length(y)
    c = t;
    t = y(:);
    t = t(:)-min(abs(t(:))); 
    zz = [ones(m,1) exp(-t/tau(2)) exp(-t*(1/tau(1)+1/p))];
    err = zz*c;
else
    t = t(:)-min(abs(t(:))); 
    y = y(:);
	t = t(isfinite(y));
    y = y(isfinite(y));
    m = length(t);
    zz = [ones(m,1) exp(-t/tau(2)) exp(-t*(1/tau(1)+1/p))];
    c = zz\y;
    z = zz*c;
    if nargin>3 && ~isempty(pic)
        if pic==1
            plot(t, y, 'ob', t, z, 'r'); drawnow
        else
            semilogy(t, y, 'o', t, z); drawnow
        end
    end

    err = sum(sum((y-z).^2./abs(z)));
    %err = sum((y-z).^2);
    %ind = y>0;
    %err = sum(y(ind).*log(y(ind)./z(ind))-y(ind)+z(ind));
end

