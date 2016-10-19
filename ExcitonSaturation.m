function [y, z] = ExcitonSaturation(x,n,kc,t)

% saturation of N excitonically coupled emitters

x = x(:);
fac = n;
y = fac*x.*(1-x).^(n-1);
for j=2:n
   fac = fac*(n-j+1)/j;
   y = y + fac*x.^j.*(1-x).^(n-j)*sum(1./(1+kc*(0:j-1)));
end
    
if nargin>3
    weight = zeros(1,n);
    weight(n) = y;
    M = zeros(n);
    M(n,n) = -y;
    if n>1
        M(n-1,n) = y;
    end
    for k=n-1:-1:2
        fac = k;
        tmp = fac*x.*(1-x).^(k-1);
        for j=2:k
            fac = fac*(k-j+1)/j;
            tmp = tmp + fac*x.^j.*(1-x).^(k-j)*sum(1./(1+kc*(0:j-1)));
        end
        M(k,k) = -tmp;
        M(k-1,k) = tmp;
        weight(k) = tmp;
    end
    M(1,1) = -x;
    weight(1) = x;
    z0 = zeros(n,1);
    z0(end) = 1;
    for j=1:length(t)
        z(j) = weight*(expm(t(j)*M)*z0);
    end
end

y = y(end);
z = z/y;

