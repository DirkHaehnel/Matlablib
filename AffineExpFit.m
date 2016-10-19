function [err, c, z] = AffineExpFit(p, td, yd, tm, ym, yx, pic)

td = td(:); tm = p(1)*tm(:); p=p(:)';
if td(1)<tm(1) 
   tm = [td(1); tm];
   ym = [ym(1,:); ym];
end
if td(end)>tm(end)
   tm = [tm; td(end)];
   ym = [ym; zeros(1,size(ym,2))];
end

zz = interp1(tm,ym,td,'cubic');

z = 0*yd;
if nargin>5 && ~isempty(yx)
    if length(p)==5
        p(5:7) = p(2:4)*p(5);
        len = 2;
    elseif length(p)==3
        p(5:7) = p(3);
        p(2:4) = p(2);
        len = 2;
    elseif length(p)==2
        p(2:4) = p(2);
        len = 1;
    else
        len = (length(p)-1)/3;
    end
    sub = 3;
else
    len = (length(p)-1);
    sub = size(yd,2);
end
for j=1:size(yd,2)
    M = [ones(size(yd,1),1) zz(:,j)];
    for k=1:len
        M = [M zz(:,j).*exp(-td/p(j+1+(k-1)*sub))];
    end
    c(:,j) = lsqnonneg(M,yd(:,j));
    z(:,j) = M*c(:,j);
end
if nargin>5 && ~isempty(yx)
    for j=1:size(yx,2)
        z(:,2+j) = [ones(size(yd,1),1) zz(:,2+j)]*sqrt(c(1:2,1).*c(1:2,2));
        M = [];
        for k=1:len
            M = [M zz(:,2+j).*exp(-td/p(4+(k-1)*3))];
        end
        cc = lsqnonneg(M,yx(:,j)-z(:,2+j));
        z(:,2+j) = z(:,2+j) + M*cc;
        yd(:,2+j) = yx(:,j);
    end
end

if nargin>6 && ~isempty(pic)
     semilogx(td,yd,td,z); drawnow
end

err = sum(sum((yd-z).^2));

