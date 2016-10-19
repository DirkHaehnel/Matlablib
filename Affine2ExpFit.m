function [err, c, z, z1, z2] = Affine2ExpFit(p, td, yd, tm, ym, yx, pic)

td = td(:); 
tm1 = p(1)*tm(:); 
tm2 = p(2)*tm(:); 
on = ym(1,:);
ze = zeros(1,size(ym,2));
if p(1)<p(2)
    tm2 = [tm1(1); tm2];
    tm1 = [tm1; tm2(end)];
    ym2 = [on; ym];
    ym1 = [ym; ze];
else
    tm1 = [tm2(1); tm1];
    tm2 = [tm2; tm1(end)];
    ym1 = [on; ym];
    ym2 = [ym; ze];
end
p=p(:)';
if td(1)<tm1(1)
   tm1 = [td(1); tm1];
   tm2 = [td(1); tm2];
   ym1 = [on; ym1];
   ym2 = [on; ym2];   
end
if td(end)>tm1(end)
   tm1 = [tm1; td(end)];
   tm2 = [tm2; td(end)];
   ym1 = [ym1; ze];
   ym2 = [ym2; ze];   
end

z1 = interp1(tm1,ym1,td,'cubic');
z2 = interp1(tm2,ym2,td,'cubic');

z = 0*yd;
if nargin>5 && ~isempty(yx)
    if length(p)==6
        p(6:8) = p(3:5)*p(6);
        len = 2;
    elseif length(p)==4
        p(6:8) = p(4);
        p(3:5) = p(3);
        len = 2;
    elseif length(p)==3
        p(3:5) = p(3);
        len = 1;
    else
        len = (length(p)-2)/3;
    end
    sub = 3;
else
    len = (length(p)-2);
    sub = 2;
end
for j=1:size(yd,2)
    M = [ones(size(yd,1),1) z1(:,j) z2(:,j)];
    for k=1:len
        M = [M exp(-td/p(j+2+(k-1)*sub))];
    end
    c(:,j) = lsqnonneg(M,yd(:,j));
    z(:,j) = M*c(:,j);
end
if nargin>5 && ~isempty(yx)
    for j=1:size(yx,2)
        z(:,2+j) = [ones(size(yd,1),1) z1(:,2+j) z2(:,2+j)]*sqrt(c(1:3,1).*c(1:3,2));
        M = [];
        for k=1:length(p)-2-2*len
            M = [M exp(-td/p(5+(k-1)*3))];
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

