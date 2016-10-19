% FCSCell

NA = 1.2;
n0 = 1.33
n
n1
lamex
over
focpos
resolution
maxnum = 1e3;

k0 = 2*pi/lamex*n0; 
k = 2*pi/lamex*n; 
k1 = 2*pi/lamex*n1; 

chimax = asin(NA/n0(1)); 
chimin = 0;

dchi = (chimax-chimin)/maxnum;
chi = (chimin:dchi:chimax)';

c0 = cos(chi);
s0 = sin(chi);
ss = n0(1)/n*s0;
cc = sqrt(1-ss.^2);
rad = over(1)*s0/sin(chimax); % relation between mode angle and entry ray lateral distance

if over(3)==0
    fac = 0*chi;    
    fac(1) = 1;
else
    if over(2)==inf
        zeta = 0;
        ww = over(3);
    else
        zeta = over(2)*lamex/pi/over(3)^2;
        ww = over(3)*sqrt(1+zeta^2);
    end
    fac = sin(chi).*sqrt(c0).*exp(-rad.^2/ww^2*(1-i*zeta));
end
    
phase1 = k*cc*zv - k0(1)*c0*focpos*ones(size(zv));
phase2 = k*cc*(2*d-zv) - k0(1)*c0*focpos*ones(size(zv));
ez1 = exp(i*phase1);
ez2 = exp(i*phase2);

