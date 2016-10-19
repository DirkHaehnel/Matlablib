function [err, c, z] = RiglerZMW(p, t, y)

% Function RiglerZMW(p, t, y) calculates the error err, amplitudes c, and fit curve z for the 
% FCS in an exponenential 1D PSF

if length(t)<length(y)
    c = t;
    t = y(:);
    p = p(:)';
    mu = sqrt(p(1)*t);
    z = [ones(length(t),1) mu/sqrt(pi)+(1-2*mu.^2)/2.*erfcx(mu) exp(-t*p(2:end))];
    err = z*c;
else
    t = t(:);
    y = y(:);
    p = p(:)';
    mu = sqrt(p(1)*t);
    z = [ones(length(t),1) mu/sqrt(pi)+(1-2*mu.^2)/2.*erfcx(mu) exp(-t*p(2:end))];
    c = lsqnonneg(z,y);
    z = z*c;
    semilogx(t,y,'--',t,z); drawnow;
    err = sum((y-z).^2./abs(z))
end
