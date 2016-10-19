function [err, c, z] = PhotophysicsSat(p,int,y,w0,lamex,a0,av,lamem)

if nargin<3 || isempty(y)
    p = p(:)';
    int = int(:);

    if length(p)==1
        if nargin>3 && ~isempty(w0)
            Nmax = 1e2;
            zmax = pi*a0^2/lamem*sqrt(-2*av^2/a0^2/log(0.99)-1); % evaluate function until amplitude function has fallen off to 1%
            z = (0.5:Nmax)'/Nmax*zmax;
            rhomax = 3;
            rho = (0.5:Nmax)/Nmax*rhomax;
            ww = LaserBeamFit(w0,lamex,z);
            kappa = D1CEF(a0,av,lamem,z);
            err = zeros(size(int));
            for j=1:length(ww)
                exc = int*w0^2/ww(j)^2*exp(-2*rho.^2);
                err = err + kappa(j)*ww(j)^2*(1-exp(-p*exc))*rho'*rhomax/Nmax;
            end
        else
            err = 1-exp(-p*int);
        end
    else
        if nargin>3 && ~isempty(w0)
            Nmax = 1e2;
            zmax = pi*a0^2/lamem*sqrt(-2*av^2/a0^2/log(0.99)-1); % evaluate function until amplitude function has fallen off to 1%
            z = (0.5:Nmax)'/Nmax*zmax;
            rhomax = 3;
            rho = (0.5:Nmax)/Nmax*rhomax;
            ww = LaserBeamFit(w0,lamex,z);
            kappa = D1CEF(a0,av,lamem,z);
            err = zeros(size(int));
            for j=1:length(ww)
                exc = int*w0^2/ww(j)^2*exp(-2*rho.^2);
                err = err + kappa(j)*ww(j)^2*(p(2)*exc./(p(1)+p(2)*exc).*(1-exp(-p(1)-p(2)*exc)))*rho'*rhomax/Nmax;
            end
        else
            err = p(2)*int./(p(1)+p(2)*int).*(1-exp(-p(1)-p(2)*int));
        end
        err(isnan(err)) = 0;
    end
else
    p = p(:)';
    int = int(:);
    y = y(:);
    
    if nargin>3 && ~isempty(w0)
        z = PhotophysicsSat(p,int,'',w0,lamex,a0,av,lamem);
    else
        z = PhotophysicsSat(p,int);
    end
    c = z\y;
    z = z*c;
    
    plot(int,y,'o',int,z)
    drawnow
    
    err = sum(abs(z-y).^2)
end


