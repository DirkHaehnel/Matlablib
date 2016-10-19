function [err, c, z] = HillPoly(p, conc, y, poly, conc0)

if length(p)>2
    conc = conc(:) - p(3);
else
    conc = conc(:);
end
if size(y,1)==1 || size(y,2)==1
    y = y(:);
end

z = conc.^p(2)./(p(1)^p(2)+conc.^p(2));
if size(z,2)>1
    for j=1:size(z,2)
        zz = [ones(size(conc)) z(:,j)];
        c(:,j) = zz\y(:,j);
        z(:,j) = zz*c(:,j);
    end
else
    zz = [ones(size(y)) z];
    if nargin>3 && ~isempty(poly)
        for k=1:poly
            zz = [zz conc.^k];
        end
    end
    
    c = zz\y;
    z = zz*c;
end

if nargin<5 || isempty(conc0)
    semilogx(conc, y, 'o', conc, z); drawnow;
    %plot(conc, y, 'o', conc, z); drawnow;
end

%err = sum(sum((y-z).^2./abs(z)));
err = sum(sum((y-z).^2));

if nargin>4 && ~isempty(conc0)
    conc0 = conc0(:);
    z = conc0.^p(2)./(p(1).^p(2)+conc0.^p(2));
    if size(c,2)>1
        for j=1:size(c,2)
            z(:,j) = [ones(size(conc0)) z(:,j)]*c(:,j);
        end
    else
        z = [ones(size(conc0)) z]*c;
    end
    semilogx(conc, y, 'o', conc0, z); drawnow;
    %plot(conc, y, 'o', conc0, z); drawnow;
end

    
