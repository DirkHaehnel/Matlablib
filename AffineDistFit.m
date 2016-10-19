function [err, c, z] = AffineDistFit(p, td, yd, tm, ym0, yx, pic)

tmp = max([0 p(1)-3*p(2)]) + (p(1) + 3*p(2) - max([0 p(1)-3*p(2)]))*(0.5:100)/100;
td = td(:);
zz = [];
if nargin>5 && ~isempty(yx)
    cnt = size(yd,2)+size(yx,2);
else
    cnt = size(yd,2);
end
for j=1:length(tmp)
    t = tm(:)*tmp(j); ym = ym0;
    if td(1)<min(t)
        t = [td(1); t];
        ym = [ym0(1,:); ym0];
    end
    if td(end)>max(t)
        t = [t; td(end)];
        ym = [ym0; zeros(1,size(ym0,2))];
    end
    if isempty(zz)
        if size(ym,2)==cnt
            zz = exp(-(tmp(j)-p(1)).^2/2/(p(2)+(p(2)==0)))*interp1(t,ym,td,'cubic');
        else
            zz = exp(-(tmp(j)-p(1)).^2/2/(p(2)+(p(2)==0)))*interp1(t,ym(:,1:cnt),td,'cubic');
        end
    else
        if size(ym,2)==cnt
            zz = zz + exp(-(tmp(j)-p(1)).^2/2/(p(2)+(p(2)==0)))*interp1(t,ym,td,'cubic');
        else
            zz = zz + exp(-(tmp(j)-p(1)).^2/2/(p(2)+(p(2)==0)))*interp1(t,ym(:,1:cnt),td,'cubic');
        end
    end
end

z = 0*yd;
for j=1:size(yd,2)
    c(:,j) = lsqnonneg([ones(size(yd,1),1) zz(:,:,j)],yd(:,j));
    z(:,j) = [ones(size(yd,1),1) zz(:,:,j)]*c(:,j);
end
if nargin>5 && ~isempty(yx)
    z(:,3) = [ones(size(yd,1),1) zz(:,:,3)]*sqrt(c(:,1).*c(:,2));
    yd(:,3) = yx;
end

if nargin>6 && ~isempty(pic)
    if pic==1
        semilogx(td,yd,td,z);
    else
        plot(td,yd,td,z); 
    end
    drawnow
end

err = sum(sum((yd-z).^2));

