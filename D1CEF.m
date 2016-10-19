function [err, c, y] = D1CEF(p,a,lam,z,amp,pic,bck)

% Function D1CEF(p,a,lam,z,amp,pic,bck) calculates the collection efficiency
% (c) Joerg Enderlein, http://www.joerg-enderlein.de (2008)

if length(p)==1
    p = [p 0];
end
z = z - p(2);
R2 = p(1).^2*(1 + (lam*z/pi/p(1)^2).^2);
y = (1-exp(-2*a^2./R2))/(1-exp(-2*a^2./p(1).^2));

if nargin>4
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