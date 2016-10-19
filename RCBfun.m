function [err,y] = RCBfun(x,c,p,k1,k2);

y = (c+k2)/2;
tmp = sqrt(y.^2 - x.*(c-x));
y = y - tmp;
err = (p-x).*(c-x-y)-k1*(x-y);

