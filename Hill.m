function [err, c, z] = Hill(p, conc, y, conc0, flag)

conc = conc(:);
if size(y,1)==1 || size(y,2)==1
    y = y(:);
end

z = conc.^p(2)./(p(1)^p(2)+conc.^p(2));
if size(z,2)>1
    for j=1:size(z,2)
        c(:,j) = [ones(size(conc)) z(:,j)]\y(:,j);
        z(:,j) = [ones(size(conc)) z(:,j)]*c(:,j);
    end
else
    %c = lsqnonneg([ones(size(y)) z],y);
    c = [ones(size(y)) z]\y;
    z = [ones(size(y)) z]*c;
end

if nargin<4 || isempty(conc0)
    if nargin<5 || isempty(flag)
        semilogx(conc, y, 'o', conc, z); drawnow;
    else
        plot(conc, y, 'o', conc, z); drawnow;
    end
end

%err = sum(sum((y-z).^2./abs(z)));
err = sum(sum((y-z).^2));

if nargin>3 && ~isempty(conc0)
    conc0 = conc0(:);
    z = conc0.^p(2)./(p(1).^p(2)+conc0.^p(2));
    if size(c,2)>1
        for j=1:size(c,2)
            z(:,j) = [ones(size(conc0)) z(:,j)]*c(:,j);
        end
    else
        z = [ones(size(conc0)) z]*c;
    end
    if nargin<5 || isempty(flag)
        semilogx(conc, y, 'o', conc0, z); drawnow;
    else
        plot(conc, y, 'o', conc0, z); drawnow;
    end
end

    
