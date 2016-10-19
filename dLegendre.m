function pm = dLegendre(n,m,x)

[pm, pm] = mLegendre(n,abs(m),x);
pm = pm.*sqrt(1-x(:).^2); 
