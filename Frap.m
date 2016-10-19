function [y, z] = frap(x,t);

x = x(:)';
k = (0.5:3e4)'/10;

tmp = -exp(-k.^2*t).*besselj(1,k);

for j=1:length(x)
    y(j) = 1+sum(tmp.*besselj(0,k*x(j)))*mean(diff(k));
end

if nargout>1
    for j=1:length(x) 
        ind = x>x(j); 
        z(j) = sum(x(ind).*y(ind)./sqrt(x(ind)-x(j))); 
    end    
end