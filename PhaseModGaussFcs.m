function [err, c, z, zz, back] = PhaseModGaussFcs(para,av,lam,delta,td,yd,expflag,bounds)

% para(1) - laser beam
% para(2) - pupil parameter

global pd
if nargin>5
    pd
end

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
ym = ym*[1 1; 1 -1];

if nargin==5
    td = td(:);
    tm = tm(:);
    if td(1)<min(tm(1))
        tm = [td(1); tm];
        ym = [ym(1,:); ym];
    end
    if td(end)>max(tm(end))
        tm = [tm; td(end)];
        ym = [ym; zeros(1,size(ym,2))];
    end
    ym = interp1(tm,ym,td,'cubic');
end

if nargin>5 && ~isempty(yd)
    td = td(:);
    if ~(size(yd,1)==length(td))
        yd = yd';
    end
    if nargin<8 || isempty(bounds)
        pdmin = zeros(size(pd));
        pdmax = [];
    else
        pdmin = bounds(1,:);
        pdmax = bounds(2,:);
    end
    pd = simplex('AffineFit',pd,pdmin,pdmax,[],[],td, yd, tm, ym, [], [], expflag);
    [err, c, z, zz, back] = AffineFit(pd, td, yd, tm, ym, [], 1, expflag);
else
    err = ym;
    c = tm;
end

