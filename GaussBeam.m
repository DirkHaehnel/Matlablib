function [err, c, x] = GaussBeamFit(p,lam,z,y,flag,pic,weight)

% fitting the Gaussian radius of an aberrated Gauss beam
if length(p)==1
    p = [p 0];
end
if length(p)>2
    x = p(1)*sqrt(1+(lam*(z-p(2))/pi/p(3)^2).^2);
else
    x = p(1)*sqrt(1+(lam*(z-p(2)+6*pi*p(3)/lam/p(1)^2)/pi/p(1)^2+924*p(3)^2/p(1)^8).^2);
end
if nargin>3 
    if nargin<5 || isempty(flag)
        c = x(:)\y(:);
    else
        c = 1;
    end
    x = x*c;
    if nargin>5 && ~isempty(pic)
        plot(z,y,'o',z,x); drawnow
    end
    if nargin>6 && ~isempty(weight)
        err = sum((y(:)-x(:)).^2./x(:).*weight(:));
    else
        err = sum((y(:)-x(:)).^2./x(:));
    end
else
    err = x;
end

