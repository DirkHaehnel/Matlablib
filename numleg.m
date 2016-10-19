function y=numleg(l,m,x)
A = legendre(l,x);
y = (-1)^m*A(m+1,1);