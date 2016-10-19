function [err, c, z] = PhaseModFcs(para,av,lam,delta,td,yd)

% para(1) - laser beam 
% para(2) - pupil parameter
% para(3) - modulation period 
% para(4) - diffusion coefficient in cm^2/s times 1e8 
% para(5:...) - exponents

para = para(:)';

Nmax = 1e2;
zmax = pi*para(2)^2/lam(2)*sqrt(-2*av^2/para(2)^2/log(0.99)-1); %evaluate function until amplitude function has fallen off to 1%
eta = (0.5:Nmax)/Nmax*zmax; %*av*pi/lam(2)*p(2);
xi = (0.5:Nmax)'/Nmax*6;
exi = exp(-xi'.^2);
t = 10.^-4*para(1)^2;
alpha = 10^0.1;
row = ones(1,Nmax);
flag = true;
j = 1;
while flag
    tm(j) = t;
    t2 = sqrt(t);
    w1 = LaserBeamFit([para(1) 0],lam(1),row'*eta-t2*xi*row).^2;
    w2 = LaserBeamFit([para(1) 0],lam(1),row'*eta+t2*xi*row).^2;
    amp = D1CEF([para(2) 0],av,lam(2),row'*eta-t2*xi*row).*D1CEF([para(2) 0],av,lam(2),row'*eta+t2*xi*row);
    ww = 8*t + w1 + w2;
    ym(j,1) = sum(exi*(amp./ww));
    ym(j,2) = sum(sum(exp(-xi.^2*row-2*delta^2./ww).*(amp./ww)));
    if ym(end,1)<1e-3*max(ym(:,1))
        flag = false;
    end
    t = alpha*t; j = j+1;
end
ym = ym/ym(1,1);

td = td(:);
tm = tm(:);
if para(4)*td(1)<min(tm(1))
    tm = [para(4)*td(1); tm];
    ym = [ym(1,:); ym];
end
if para(4)*td(end)>max(tm(end))
    tm = [tm; para(4)*td(end)];
    ym = [ym; zeros(1,size(ym,2))];
end
cd = cos(2*pi*td/para(3));
ym = [interp1(tm/para(4),ym(:,1),td).*(1 + 0.5*cd), ...
    interp1(tm/para(4),ym(:,2),td).*(1 - 0.5*cd) exp(-td*(1./para(5:end)))];

c = lsqnonneg([ones(length(td),1) ym],yd);
z = [ones(length(td),1) ym]*c;
semilogx(td,yd,td,z); drawnow
err = sum((z-yd).^2./z./td)

