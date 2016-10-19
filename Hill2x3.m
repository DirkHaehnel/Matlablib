function [err, c2, c3, z2, z3] = Hill2x3(p, conc, y2, y3, conc0)

conc = conc(:);
if size(y2,1)==1 || size(y2,2)==1
    y2 = y2(:);
end
if size(y3,1)==1 || size(y3,2)==1
    y3 = y3(:);
end

z2 = conc.^p(2)./(p(1)^p(2)*p(3)^2+conc.^p(2));
z3 = conc.^p(2)./(p(1)^p(2)*p(3)^3+conc.^p(2));
if size(z2,2)>1
    for j=1:size(z2,2)
        c2(:,j) = [ones(size(conc)) z2(:,j)]\y2(:,j);
        z2(:,j) = [ones(size(conc)) z2(:,j)]*c2(:,j);
        c3(:,j) = [ones(size(conc)) z3(:,j)]\y3(:,j);
        z3(:,j) = [ones(size(conc)) z3(:,j)]*c3(:,j);
    end
else
%     c2 = lsqnonneg([ones(size(y2)) z2],y2);
%     z2 = [ones(size(y2)) z2]*c2;
%     c3 = lsqnonneg([ones(size(y3)) z3],y3);
%     z3 = [ones(size(y3)) z3]*c3;
    c2 = lsqnonneg(z2,y2);
    z2 = z2*c2;
    c3 = lsqnonneg(z3,y3);
    z3 = z3*c3;
end

if nargin<5 || isempty(conc0)
    semilogx(conc, y2, 'o', conc, y3, 'v', conc, z2, conc, z3); drawnow;
end

%err = sum(sum((y-z).^2./abs(z)));
err = sum(sum((y2-z2).^2+(y3-z3).^2));

if nargin>4 && ~isempty(conc0)
    conc0 = conc0(:);
    z2 = conc0.^p(2)./(p(1).^p(2)*p(3)^2+conc0.^p(2));
    z3 = conc0.^p(2)./(p(1).^p(2)*p(3)^3+conc0.^p(2));
    if size(c2,2)>1
        for j=1:size(c2,2)
            z2(:,j) = [ones(size(conc0)) z2(:,j)]*c2(:,j);
            z3(:,j) = [ones(size(conc0)) z3(:,j)]*c3(:,j);
        end
    else
        z2 = [ones(size(conc0)) z2]*c2;
        z3 = [ones(size(conc0)) z3]*c3;
    end
    semilogx(conc, y2, 'o', conc, y3, 'x', conc0, z2, conc0, z3); drawnow;
end

    
