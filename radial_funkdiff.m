function radfunkdiff = radial_funkdiff(n,l,x)

%calculates d/dr of R_nl(r)
%uses radial_funk.m, lagr.m , norm_fac.m
%assumes z=1 (charge) and a =1 (lenght unit).
% We have included a factor x.*x. from the radial integration
% to prevent difficulties i origo.
z=1;a=1;
x1=((2*z)/(n*a)).*x;
%term1and2 = l*radial_funk(n,l,x)/x - radial_funk(n,l,x)/n; 
%term3= -norm_fac(n,l).*exp(-x1/2).*x1.^l.*lagr_diff(n,l,x1)*2/n;
term1and2 = (l - (x/n)).*x.*radial_funk(n,l,x);
%term1and2 = l*x.*radial_funk(n,l,x) - x.*x.*radial_funk(n,l,x)/n; 
term3= -x.*x.*norm_fac(n,l).*exp(-x1/2).*x1.^l.*lagr_diff(n,l,x1)*2/n;
radfunkdiff = term1and2 + term3;
 
