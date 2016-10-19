function [err, c, z] = Gauss2Fcs(para,av,lam,delta,td,yd,yx,expflag,multiflag)

% p(1) - laser beam 
% p(2) - detection

global pd
if isempty(pd)
    disp('please define global variable pd');
    return
end

Nmax = 1e2;
zmax = pi*para(3)^2/lam(2)*sqrt(-2*av^2/para(2)^2/log(0.99)-1); %evaluate functionuntil amplitude function has fallen off to 1%
eta = (0.5:Nmax)/Nmax*zmax; %*av*pi/lam(2)*p(2);
xi = (0.5:Nmax)'/Nmax*6;
tm = 10.^(-4:0.1:3)*para(1)^2;
row = ones(1,Nmax);
ym = zeros(size(tm));
for j=1:length(tm)
    t = tm(j);
    t2 = sqrt(t);
    w1 = LaserBeamFit([para(1) 0 para(2)],lam(1),row'*eta-t2*xi*row).^2;
    w2 = LaserBeamFit([para(1) 0 para(2)],lam(1),row'*eta+t2*xi*row).^2;
    if length(para)==3
        para(4) = para(3);
    end
    amp = D1CEF([para(3) 0 para(4)],av,lam(2),row'*eta-t2*xi*row).*D1CEF([para(3) 0 para(4)],av,lam(2),row'*eta+t2*xi*row);
    ww = 8*t + w1 + w2;
    ym(j,1) = sum(sum(exp(-xi'.^2)*(amp./ww)));
    if nargin>3 && ~isempty(delta)
        ym(j,2) = sum(sum(exp(-xi.^2*row-2*delta^2./ww).*(amp./ww)));
    end
end
ym = ym/ym(1,1);

if nargin>5
    td = td(:);
    if ~(size(yd,1)==length(td))
        yd = yd';
    end
    if nargin<8 || isempty(expflag)
        if nargin<7 || isempty(yx)
            if nargin<9 || isempty(multiflag)
                pd = simplex('AffineFit',pd,0,[],[],[],td, yd, tm, ym);
                [err, c, z] = AffineFit(pd, td, yd, tm, ym, [], 1);
            else
                if any(size(multiflag)>1)
                    pdmin = multiflag(1,:);
                    pdmax = multiflag(2,:);
                else
                    pdmin = zeros(size(pd));
                    pdmax = [];
                end
                pd = simplex('AffineFit',pd,pdmin,pdmax,[],[],td, yd, tm, ym);
                [err, c, z] = AffineFit(pd, td, yd, tm, ym, [], 1);
            end
        else
            pd = simplex('AffineFit',1,0,[],[],[],td, yd, tm, ym(:,[1 1 2]), yx);
            [err, c, z] = AffineFit(pd, td, yd, tm, ym(:,[1 1 2]), yx, 1);
        end
    else
        if nargin<7 || isempty(yx)
            if nargin<9 || isempty(multiflag)
                pd = simplex('AffineExpFit',pd,zeros(size(pd)),[],[],[],td, yd, tm, ym);
                [err, c, z] = AffineExpFit(pd, td, yd, tm, ym, [], 1);
            else
                if any(size(multiflag)>1)
                    pdmin = multiflag(1,:);
                    pdmax = multiflag(2,:);
                else
                    pdmin = zeros(size(pd));
                    pdmax = [];
                end
                pd = simplex('Affine2ExpFit',pd,pdmin,pdmax,[],[],td, yd, tm, ym);
                [err, c, z] = Affine2ExpFit(pd, td, yd, tm, ym, [], 1);
            end
        else
            if nargin<9 || isempty(multiflag)
                pd = simplex('AffineExpFit',pd,zeros(size(pd)),[],[],[],td, yd, tm, ym(:,[1 1 2]), yx);
                [err, c, z] = AffineExpFit(pd, td, yd, tm, ym(:,[1 1 2]), yx, 1);
            else
                if all(size(multiflag)>1)
                    pdmin = multiflag(1,:);
                    pdmax = multiflag(2,:);
                else
                    pdmin = zeros(size(pd));
                    pdmax = [];
                end
                pd = simplex('Affine2ExpFit',pd,pdmin,pdmax,[],[],td, yd, tm, ym(:,[1 1 2]), yx);
                [err, c, z] = Affine2ExpFit(pd, td, yd, tm, ym(:,[1 1 2]), yx, 1);
            end
        end
    end
else
    err = ym;
    c = tm;
end

