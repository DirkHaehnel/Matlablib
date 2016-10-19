function [err, c, z, zz, back] = Gauss3Fcs(para,av,lam,td,yd,expflag,bounds)

% para(1) - laser beam 
% para(2) - pupil parameter

if nargin>4
    global pd
    pd
end

if 0
    Nxmax = 2e2;
    Nzmax = 1e3;
    zmax = pi*para(2)^2/lam(2)*sqrt(-2*av^2/para(2)^2/log(0.9*(1-exp(-2*av^2./para(2).^2)))-1); %evaluate function until amplitude function has fallen off to 1%
    eta = (0.5:Nzmax)'/Nzmax*zmax; %*av*pi/lam(2)*p(2);
    xi = (0.5:Nxmax)/Nxmax*5;
    xi2 = xi.^2;
    t = 1e-6*para(1)^2;
    alpha = 10^0.1;
    row = ones(1,length(xi));
    col = ones(length(eta),1);
    flag = true;
    j = 1;
    w = LaserBeamFit([para(1) 0],lam(1),eta);
    w2 = w.^2;
    amp = D1CEF([para(2) 0],av,lam(2),eta);
    exi = ((amp./w2)*row).*exp(-2*w2*xi2);
    while flag
        tm(j) = t;
        ww = 8*t + w2;
        ff = ((amp./ww)*row).*exp(-2*(w2./ww)*xi2);
        ff = (exp(-(col*eta'-eta*col').^2/4/t)+exp(-(col*eta'+eta*col').^2/4/t))*ff;
        %ff = mconv2([flipud(ff); ff],exp(-[-flipud(eta); eta].^2/4/t));
        %ff = ff(1:end/2,:);
        ym(j,1) = sum(sum((w2*xi).*(exi.*ff.^2)))/t;
        if ym(end,1)<1e-5*max(ym(:,1))
            flag = false;
        end
        t = alpha*t; j = j+1;
    end
else
    t = 1e-10*para(1)^2;
    alpha = 10^0.1;
    flag = true;
    j = 1;
    while flag
        tm(j) = t;
        ym(j,1) = 1./(64*t^2+32*para(1)^2*t+3*para(1)^4)/sqrt(64*t^2+32*para(2)^2*t+3*para(2)^4);
        if ym(end,1)<1e-5*max(ym(:,1))
            flag = false;
        end
        t = alpha*t; j = j+1;
    end    
end
ym = ym/ym(1,1);


if nargin==4
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
    tm = td;
end

if nargin>4 && ~isempty(yd)
    td = td(:);
    if ~(size(yd,1)==length(td))
        yd = yd';
    end
    if nargin<7 || isempty(bounds)
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

