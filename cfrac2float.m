function x = cfrac2float(b,a);

eps = 1e-30;

if b(1)>0
    f(1) = b(1);
else
    f(1) = eps;
end

c(1) = f(1);
d(1) = 0;

for j=1:length(a)
    d(j+1) = b(j+1) + a(j)*d(j);
    if d(j+1)==0 d(j+1)=eps; end
    c(j+1) = b(j+1) + a(j)/c(j);
    if c(j+1)==0 c(j+1)=eps; end
    d(j+1) = 1/d(j+1);
    delta(j) = c(j+1)*d(j+1);
    f(j+1) = f(j)*delta(j);
end

x = f(end);


