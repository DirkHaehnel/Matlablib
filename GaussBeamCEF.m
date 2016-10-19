function [err, c, y] = GaussBeamCEF(p,w,a,lam,z,amp,pic,bck)

% Function D1CEF(p,a,lam,z,amp,pic,bck) calculates the collection efficiency

if length(p)==1
    p = [p 0];
end
R2 = p(1).^2*(1 + (lam(2)*(z-p(2))/pi/p(1)^2).^2);
R0 = p(1).^2*(1 + (lam(2)/pi/p(1)^2).^2);
y = (1-exp(-2*a^2./R2))/(1-exp(-2*a^2./R0));
if length(p)==3
    nmax = 10;
    xi = (w^2+i*lam(1)*z/pi)./8/sqrt(i*p(3));
    tmp = 1;
    for j=1:length(nmax)
        tmp = tmp + (-1)^j*gamma(2*j+1)/gamma(j+1)./(2*xi).^(2*j);
    end;
    y = y.*abs(tmp).^2;
end

if nargin>5
    z = z(:); amp = amp(:); y = y(:);
    if nargin>6 && ~isempty(bck)
    	c = [ones(size(z)) y]\amp;
        y = [ones(size(z)) y]*c;
    else
        c = y\amp;
        y = y*c;
    end

    if nargin>5 && ~isempty(pic)
        plot(z,amp,'o',z,y); drawnow
    end

    err = sum((amp-y).^2./y);
else
    err = y;
end