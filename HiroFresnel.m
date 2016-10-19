function [Rp, Rs, Tp, Ts] = HiroFresnel(theta, n1, n2)

if (theta==0.)
    Rp = (-n1+n2)./(n1+n2);
    Rs =  (n1-n2)./(n1+n2);
    Tp = 2.*n1./(n1+n2);
    Ts = 2.*n1./(n1+n2);
elseif (theta == pi/2)
    Rp = 1;
    Rs = 1;
    Tp = 0;
    Ts = 0;
else
    phi1 = theta;
    phi2 = asin(n1.*sin(phi1)./n2);
    Rp =  tan(phi1 - phi2)./tan(phi1 + phi2);
    Rs = -sin(phi1 - phi2)./sin(phi1 + phi2);
    Tp = 2.*cos(phi1).*sin(phi2)./sin(phi1 + phi2)./cos(phi1 - phi2);
    Ts = 2.*cos(phi1).*sin(phi2)./sin(phi1 + phi2);
end

