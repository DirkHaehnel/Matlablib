function [err, c, x] = GaussBeamFit(p,lam,z,y,pic,weight)

% fitting the Gaussian radius of an aberrated Gauss beam
if length(p)==1
    p = [p 0];
end
z = z - p(2);
if length(p)==2 || p(3)==0 
    x = p(1)*sqrt(1+(lam*z/pi/p(1)^2).^2);
else
    x = p(1)*sqrt(1+(lam*(z+6*pi*p(3)/lam/p(1)^2)/pi/p(1)^2).^2+924*p(3)^2/p(1)^8);
    nmax = 5;
    xi = (p(1)^2+i*lam*z/pi)./8/sqrt(i*p(3));
    %     amp = 1;
    %     for j=1:length(nmax)
    %         amp = amp + (-1)^j*gamma(2*j+1)/gamma(j+1)./(2*xi).^(2*j);
    %     end;
    %     amp = abs(amp).^2;
    tmp = p(1)*sqrt(1+(lam*z/pi/p(1)^2).^2); 
    amp = pi*p(1)^2/64/abs(p(3))*exp(p(1)^2*lam*z/16/pi/p(3)).*abs(1-erfz(xi)).^2.*x.^2;
end
if nargin>3 
    if nargin>4 && ~isempty(pic)
        plot(z,y,'o',z,x); drawnow
    end
    if nargin>5 && ~isempty(weight)
        err = sum((y(:)-x(:)).^2./x(:).*weight(:));
    else
        err = sum((y(:)-x(:)).^2./x(:));
    end
else
    err = x;
    if exist('amp')==1
        c = amp;
    end
end

