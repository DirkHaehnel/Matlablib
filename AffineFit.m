function [err, c, z, zz, back] = AffineFit(p, td, yd, tm, ym, yx, pic, expflag, shiftflag, weight, fitfun)

% Function AffineFit for affine fitting of a model against data
% (c) Joerg Enderlein, http://www.joerg-enderlein.de (2010)

p = p(:); 
td = td(:);
if nargin>8 && ~isempty(shiftflag)
    tm = tm(:)-p(1);
    p(1) = [];
else
    tm = tm(:);
end
if nargin<10 || isempty(weight)
    weight = 1;
end
% define whether to use least-square or nonnegatve least-square fitting:
if nargin<11 || isempty(fitfun)
    fitfun = 'lsqnonneg';
else
    fitfun = 'lscov';
end
if size(yd,2)>size(yd,1)
    yd=yd.';
end
if size(ym,2)>size(ym,1)
    ym=ym.';
end
if nargin>5 && ~isempty(yx)
    if size(yx,2)>size(yx,1)
        yx=yx.';
    end
end

if ~(nargin<8 || isempty(expflag))
    if expflag(1)>0
        tau = 1./p(end-expflag(1)+1:end); tau = tau(:)';
        p = p(1:end-expflag(1));
    end
    if length(expflag)==2 && expflag(2)>0
        ww = p(end-expflag(2)+1:end);
        p = p(end-2*expflag(2)+1:end-expflag(2));
        w = (-3:0.3:3)'; tmp = [];
        for j=1:expflag(2)
            tmp = [tmp p(j)*exp(w*ww(j))];
        end
        p = [tmp(:); p(1:end-2*expflag(2))];
    else
        expflag(2) = 0;
    end
else
    expflag = [0 0];
end

td = td(:); t = tm(:)*abs(p(:)'); 
if td(1)<max(t(1,:)) 
   t = [td(1)*ones(1,length(p)); t];
   ym = [ym(1,:); ym];
end
if td(end)>min(t(end,:))
   t = [t; td(end)*ones(1,length(p))];
   ym = [ym; ym(end,:)];
end

zz = [];
if nargin>5 && ~isempty(yx)
    cnt = size(yd,2)+size(yx,2);
else
    cnt = size(yd,2);
end
if size(ym,2)==cnt
    for k=1:length(p)
        zz(:,k,:) = interp1(t(:,k),ym,td,'cubic');
    end
else
    for k=1:length(p)
        [tmp, ind] = unique(t(:,k));
        size(zz)
        size(ym)
        [(k-1)*cnt+1 k*cnt] 
        zz(:,k,:) = interp1(tmp,ym(ind,(k-1)*cnt+1:k*cnt),td,'cubic');
    end
end

if expflag(2)>0
    w = exp(-w.^2/2);
    for j=1:expflag(2)
        zz(:,j,:) = w(1)*zz(:,j,:);
        for k=2:length(w)
            zz(:,j,:) = zz(:,j,:) + w(k)*zz(:,j+k-1,:);
        end
        zz(:,j+1:j+length(w)-1,:) = [];
    end
end

if expflag(1)<0 % no flat background
    back = []; %zeros(size(yd,1),1);
elseif expflag(1)==0 % only flat background
    back = ones(size(yd,1),1); 
else % flat background + exponent(s)
    back = [ones(size(yd,1),1) exp(-td*tau)];
end

z = 0*yd;
for j=1:size(yd,2)
    eval(['c(:,j) = ' fitfun '([back zz(:,:,j)],yd(:,j));']);
end
% make sure that relative concentraion of several diffusing components
% within two foci is the same:
if size(zz,2)>1
    tmp = mean(c(end-size(zz,2)+1:end,:)./(ones(size(zz,2),1)*sum(c(end-size(zz,2)+1:end,:))),2);
    c(end-size(zz,2)+1:end,:) = tmp*sum(c(end-size(zz,2)+1:end,:),1);
end
for j=1:size(yd,2)
    z(:,j) = [back zz(:,:,j)]*c(:,j);
end
if nargin>5 && ~isempty(yx)
    for j=1:size(yx,2)
        eval(['c(:,2+j) = ' fitfun '([back zz(:,:,2+j)],yx(:,j));']);
    end
    if size(zz,2)>1
        tmp = mean(c(end-size(zz,2)+1:end,:)./(ones(size(zz,2),1)*sum(c(end-size(zz,2)+1:end,:))),2);
        c(end-size(zz,2)+1:end,:) = tmp*sum(c(end-size(zz,2)+1:end,:),1);
    end
    if size(yx,2)==4
        tmp = sqrt(c(end-size(zz,2)+1:end,1).*c(end-size(zz,2)+1:end,2)./c(end-size(zz,2)+1:end,3)./c(end-size(zz,2)+1:end,5));
        c(end-size(zz,2)+1:end,3) = c(end-size(zz,2)+1:end,3).*tmp/2;
        c(end-size(zz,2)+1:end,5) = c(end-size(zz,2)+1:end,5).*tmp/2;
        tmp = sqrt(c(end-size(zz,2)+1:end,1).*c(end-size(zz,2)+1:end,2)./c(end-size(zz,2)+1:end,4)./c(end-size(zz,2)+1:end,6));
        c(end-size(zz,2)+1:end,4) = c(end-size(zz,2)+1:end,4).*tmp/2;
        c(end-size(zz,2)+1:end,6) = c(end-size(zz,2)+1:end,6).*tmp/2;
        c(isnan(c)) = 0;
    end
    for j=1:size(yx,2)
        z(:,2+j) = [back zz(:,:,2+j)]*c(:,2+j);
        yd(:,2+j) = yx(:,j);
    end
    if size(yx,2)==4
        z(:,3) = sum(z(:,[3 5]),2);
        z(:,4) = sum(z(:,[4 6]),2);
        z(:,5:6) = [];
        yd(:,3) = sum(yd(:,[3 5]),2);
        yd(:,4) = sum(yd(:,[4 6]),2);
        yd(:,5:6) = [];
        if size(weight,1)==size(yd,1)
            weight(:,3) = sum(weight(:,[3 5]),2);
            weight(:,4) = sum(weight(:,[4 6]),2);
            weight(:,5:6) = [];
        end
    end
end

if nargin>6 && ~isempty(pic)
    if pic==1
        h = semilogx(td,yd,'o'); hold on; semilogx(td,z); hold off
    else
        plot(td,yd,'o'); hold on; plot(td,z); hold off
    end
    drawnow
end

err = sum(sum(weight.*(yd-z).^2./abs(z)));

