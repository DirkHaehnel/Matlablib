function A=norm_fac(n,l)
% A=norm_fac(n,l)
% calculate normalication factor for radial wave function
a=1;
Z=1;

A=sqrt((((2*Z)/(n*a))^3*prod(1:(n-l-1)))/(2*n*(prod(1:(n+l))^3)));
