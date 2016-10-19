function [err, c, z] = Sigmoid(p, conc, y, conc0)

conc = conc(:);
y = y(:);
if length(p)<3
    p(3) = 1;
end

z = 1./(1+p(3)*exp(-p(2)*(conc-p(1))));
c = [ones(size(y)) z]\y;
z = [ones(size(y)) z]*c;

plot(conc, y, 'o', conc, z); drawnow;

err = abs(sum((y-z).^2./abs(z)));

if nargin>3 && ~isempty(conc0)
    conc0 = conc0(:);
    z = 1./(1+p(3)*exp(-p(2)*(conc0-p(1))));
    z = [ones(size(conc0)) z]*c;
    plot(conc, y, 'o', conc0, z); drawnow;
end

    
