function z = RotoDiffEmExPlot(t,gamma,taurot,ww)

drot = 1/6/taurot;
t = t(:);

if nargin>3 && ~isempty(ww)
    load ConfoAnisoCoefficients
    z = zeros(length(t),3);
    for L=0:4
        tmpo = 0; tmop = 0; tmpp = 0;
        for M=0:L
            tmpo = tmpo + (1+sign(M))*interp2(0:15:90, w0, squeeze(permute(po(L+1,M+1,:,:),[3 4 1 2])), gamma, ww, 'cubic');
            tmop = tmop + (1+sign(M))*interp2(0:15:90, w0, squeeze(permute(op(L+1,M+1,:,:),[3 4 1 2])), gamma, ww, 'cubic');
            tmpp = tmpp + (1+sign(M))*interp2(0:15:90, w0, squeeze(permute(pp(L+1,M+1,:,:),[3 4 1 2])), gamma, ww, 'cubic');
        end
        z = z + [tmpo*exp(-L*(L+1)*drot*t), tmpp*exp(-L*(L+1)*drot*t), tmop*exp(-L*(L+1)*drot*t)];
    end
else
    z(:,1) = (4.*(2 + cos(2.*gamma))^2)./225 + (16.*(1 + 3.*cos(2.*gamma))^2)./(11025.*exp(20.*drot.*t)) + (13 + 11.*cos(2.*gamma))^2./(2205.*exp(6.*drot.*t)) + (16.*sin(gamma)^4)./(2205.* exp(20.*drot.*t)) + (4.*sin(gamma)^4)./(735.*exp(6.*drot.*t)) + (32.*sin(2.*gamma)^2)./(2205.*exp(20.*drot.*t)) +  (12.*sin(2.*gamma)^2)./(245.*exp(6.*drot.*t));
    z(:,2) = (-2.*(-3 + cos(2.*gamma)).*(2 + cos(2.*gamma)))./225 - (8.*(1 + 3.*cos(2.*gamma))^2)./(11025.*exp(20.*drot.*t)) + (269 + 44.*cos(2.*gamma) - 121.*cos(4.*gamma))./(8820.*exp(6.*drot.*t)) - (8.*sin(gamma)^4)./(2205.*exp(20.*drot.*t)) - (2.*sin(gamma)^4)./(735.*exp(6.*drot.*t)) - (16.*sin(2.*gamma)^2)./(2205.*exp(20.*drot.*t)) - (6.*sin(2.*gamma)^2)./(245.*exp(6.*drot.*t));
    z(:,3) = (-3 + cos(2.*gamma))^2./225 + (1 + 3.*cos(2.*gamma))^2./(1225.*exp(20.*drot.*t)) + (151 - 156.*cos(2.*gamma) + 37.*cos(4.*gamma))./(4410.*exp(6.*drot.*t)) + sin(gamma)^4./(245.*exp(20.*drot.*t)) - sin(gamma)^4./(63.*exp(12.*drot.*t)) + (124.*sin(gamma)^4)./(2205.*exp(6.*drot.*t)) + (2.*sin(2.*gamma)^2)./(245.*exp(20.*drot.*t)) - (2.*sin(2.*gamma)^2)./(315.*exp(12.*drot.*t)) + (52.*sin(2.*gamma)^2)./(2205.*exp(6.*drot.*t));
end

plot(t,(z./z(end,1))*diag([1 3 9]));