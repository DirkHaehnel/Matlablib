function [err, c, z, zz, back] = GaussFcs(para,av,lam,delta,td,yd,yx,expflag,bounds,flag2d)

% Function GaussFcs used by 2fFCS for fitting 2fFCS data
% (c) Joerg Enderlein, http://www.joerg-enderlein.de (2008)

% para(1) - laser beam 
% para(2) - pupil parameter
% para(3:5) - velocity vector as {vx, vy, vz}

global pd
if nargin>5
    pd
end
if nargin<10
    flag2d = 0;
end

if flag2d
    t = 10.^-4*para(1)^2;
    alpha = 10^0.1;
    flag = true;
    j = 1;
    while flag
        tm(j) = t;
        w1 = 4*t + para(1).^2;
        w2 = 4*t + para(2).^2;
        ym(j,1) = 1./w1;
        ym(j,2) = 1./w2;
        ym(j,3) = 2*exp(-2*delta^2./(w1+w2))./(w1+w2);
        if ym(end,1)<1e-3*max(ym(:,1))
            flag = false;
        end
        t = alpha*t; j = j+1;
    end
    ym = ym/max(ym(1,1:2));
else
    switch length(para)
        case {3,4,5}
            if length(para)==3
                para(4:5) = 0;
            elseif length(para)==4
                para(5) = 0;
            end
            para(3:4) = para(3:4)/1e3; % conversion from mum/mus to nm/mus
            para(5) = para(5)*sqrt(pd(1)/4e-8)/1e3;
            Nmax = 1e2;
            zmax = pi*para(2)^2/lam(2)*sqrt(-2*av^2/para(2)^2/log(0.99)-1); %evaluate function until amplitude has fallen off to 1%
            eta = (0.5:Nmax)/Nmax*zmax; %*av*pi/lam(2)*p(2);
            xi = (0.5:Nmax)'/Nmax*6;
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
                ym(j,1) = sum(sum(exp(-(xi-para(5)*t2).^2*row-2*sum((para(3:4)*t).^2)./ww).*(amp./ww)));
                
                %                 F = @(y) exp(-(xi-para(5)*t2).^2).*D1CEF([para(2) 0],av,lam(2),y-t2*xi).*...
                %                     D1CEF([para(2) 0],av,lam(2),y+t2*xi)./(8*t + LaserBeamFit([para(1) 0],lam(1),y-t2*xi).^2 + ...
                %                     LaserBeamFit([para(1) 0],lam(1),y+t2*xi).^2).*exp(-2*((para(3)*t)^2+(para(4)*t).^2)./(8*t + ...
                %                     LaserBeamFit([para(1) 0],lam(1),y-t2*xi).^2 + LaserBeamFit([para(1) 0],lam(1),y+t2*xi).^2));
                %                 ym(j,1) = sum(quadv(F,0,zmax));

                if nargin>3 && ~isempty(delta)
                    ym(j,2) = ym(j,1);

                    ym(j,3) = sum(sum(exp(-(xi-para(5)*t2).^2*row-2*((delta-para(3)*t)^2+(para(4)*t)^2)./ww).*(amp./ww)));
                    ym(j,4) = sum(sum(exp(-(xi-para(5)*t2).^2*row-2*((delta+para(3)*t)^2+(para(4)*t)^2)./ww).*(amp./ww)));

                    %                     F = @(y) exp(-(xi-para(5)*t2).^2).*D1CEF([para(2) 0],av,lam(2),y-t2*xi).*...
                    %                         D1CEF([para(2) 0],av,lam(2),y+t2*xi)./(8*t + LaserBeamFit([para(1) 0],lam(1),y-t2*xi).^2 + ...
                    %                         LaserBeamFit([para(1) 0],lam(1),y+t2*xi).^2).*exp(-2*((delta-para(3)*t)^2+(para(4)*t).^2)./(8*t + ...
                    %                         LaserBeamFit([para(1) 0],lam(1),y-t2*xi).^2 + LaserBeamFit([para(1) 0],lam(1),y+t2*xi).^2));
                    %                     ym(j,3) = sum(quadv(F,0,zmax,1e-10));
                end
                if ym(end,1)<1e-3*max(ym(:,1))
                    flag = false;
                end
                t = alpha*t; j = j+1;
            end
            ym = ym/ym(1,1);
        case 2
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
                if nargin>3 && ~isempty(delta)
                    ym(j,2) = ym(j,1);
                    ym(j,3) = sum(sum(exp(-xi.^2*row-2*delta^2./ww).*(amp./ww)));
                end
                if ym(end,1)<1e-3*max(ym(:,1))
                    flag = false;
                end
                t = alpha*t; j = j+1;
            end
            ym = ym/ym(1,1);
        case 1
            t = 10.^-4*para(1)^2;
            alpha = 10^0.1;
            flag = true;
            j = 1;
            while flag
                tm(j) = t;
                ww = 4*t + para(1).^2;
                ym(j,1) = 1./ww;
                if nargin>3 && ~isempty(delta)
                    ym(j,2) = ym(j,1);
                    ym(j,3) = exp(-delta^2./ww)./ww;
                end
                if ym(end,1)<1e-3*max(ym(:,1))
                    flag = false;
                end
                t = alpha*t; j = j+1;
            end
            ym = ym/ym(1,1);
    end
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
end

if nargin>5 && ~isempty(yd)
    td = td(:);
    if ~(size(yd,1)==length(td))
        yd = yd';
    end
    if nargin<9 || isempty(bounds)
        pdmin = zeros(size(pd));
        pdmax = [];
    else
        pdmin = bounds(1,:);
        pdmax = bounds(2,:);
    end
    if nargin<7
        yx = [];
    end
    pd = simplex('AffineFit',pd,pdmin,pdmax,[],[],td, yd, tm, ym, yx, [], expflag);
    [err, c, z, zz, back] = AffineFit(pd, td, yd, tm, ym, yx, 1, expflag);
else
    err = ym;
    c = tm;
end

