function radfunc=radial_funk(n,l,x)

%radfunc=radial_funk(n,l,x)
%Beregn radiell funksjoner for Hydrogenets boelgefunksjoner.
%Bruk: lagr.m , norm_fac.m
z=1;a=1;
%dette bestemmer ladning og enhet.
x1=((2*z)/(n*a)).*x;
radfunc=-norm_fac(n,l)*exp(-x1/2).*x1.^l.*lagr(n,l,x1);
