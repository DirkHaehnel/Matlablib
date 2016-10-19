function [err, c, z] = GaussFlowFcs(para,av,lam,delta,td,yd,yx,expflag,multiflag)

% p(1) - laser beam 
% p(2) - pupil parameter
% p(3:5) - velocity

if nargin>5
    global pd
    pd
end

para = [para(:); zeros(5-length(para),1)];

Nmax = 1e2;
zmax = pi*para(2)^2/lam(2)*sqrt(-2*av^2/para(2)^2/log(0.99)-1); %evaluate function until amplitude function has fallen off to 1%
eta = (0.5:Nmax)/Nmax*zmax; %*av*pi/lam(2)*p(2);
xi = (0.5:Nmax)'/Nmax*6;
t = 10.^-4*para(1)^2;
alpha = 10^0.05;
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
    ym(j,1) = sum(sum(exp(-(xi-para(5)*t2/2).^2*row-2*t.^2*(para(3).^2+para(4).^2)./ww).*(amp./ww)));
    if nargin>3 && ~isempty(delta)
        ym(j,2) = ym(j,1);
        ym(j,3) = sum(sum(exp(-(xi-para(5)*t2/2).^2*row-2*((delta-para(3)*t)^2+para(4)^2*t^2)./ww).*(amp./ww)));
    end
    if ym(end,1)<1e-3*max(ym(:,1))
        flag = false;
    end
    t = alpha*t; j = j+1;
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
                pd = simplex('AffineFit',pd,0*pd,[],[],[],td, yd, tm, ym);
            else
                if any(size(multiflag)>1)
                    pdmin = multiflag(1,:);
                    pdmax = multiflag(2,:);
                else
                    pdmin = zeros(size(pd));
                    pdmax = [];
                end
                pd = simplex('AffineFit',pd,pdmin,pdmax,[],[],td, yd, tm, ym);
            end
            [err, c, z] = AffineFit(pd, td, yd, tm, ym, [], 1);
        else
            if nargin<9 || isempty(multiflag)
                pd = simplex('AffineFit',pd,0*pd,[],[],[],td, yd, tm, ym, yx);
            else
                if any(size(multiflag)>1)
                    pdmin = multiflag(1,:);
                    pdmax = multiflag(2,:);
                else
                    pdmin = zeros(size(pd));
                    pdmax = [];
                end
                pd = simplex('AffineFit',pd,pdmin,pdmax,[],[],td, yd, tm, ym, yx);
            end
            [err, c, z] = AffineFit(pd, td, yd, tm, ym, yx, 1);
        end
    else
        if nargin<7 || isempty(yx)
            if nargin<9 || isempty(multiflag)
                pd = simplex('AffineExpFit',pd,0*pd,[],[],[],td, yd, tm, ym);
                [err, c, z] = AffineExpFit(pd, td, yd, tm, ym, [], 1);
            else
                if any(size(multiflag)>1)
                    pdmin = multiflag(1,:);
                    pdmax = multiflag(2,:);
                else
                    pdmin = zeros(size(pd));
                    pdmax = [];
                end
                pd = simplex('AffineExpFit',pd,pdmin,pdmax,[],[],td, yd, tm, ym);
                [err, c, z] = AffineExpFit(pd, td, yd, tm, ym, [], 1);
            end
        else
            if nargin<9 || isempty(multiflag)
                pd = simplex('AffineExpFit',pd,0*pd,[],[],[],td, yd, tm, ym, yx);
                [err, c, z] = AffineExpFit(pd, td, yd, tm, ym, yx, 1);
            else
                if all(size(multiflag)>1)
                    pdmin = multiflag(1,:);
                    pdmax = multiflag(2,:);
                else
                    pdmin = zeros(size(pd));
                    pdmax = [];
                end
                pd = simplex('AffineExpFit',pd,pdmin,pdmax,[],[],td, yd, tm, ym, yx);
                [err, c, z] = AffineExpFit(pd, td, yd, tm, ym, yx, 1);
            end
        end
    end
else
    err = ym;
    c = tm;
end

