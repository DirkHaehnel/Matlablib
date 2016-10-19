function [err, y, c, z] = Hat(p, x, ind);

[t1, t2] = meshgrid(1:size(x,1),1:size(x,2));

z = exp(-((t1-p(1)).^2+(t2-p(2)).^2)/p(3)^2/2);
x = x(logical(ind));
z = z(logical(ind));
len = prod(size(x));

z = [ones(len,1) reshape(z,len,1)];

c = z\x;
y = z*c;

err = sum((x-y).^2./abs(y))

