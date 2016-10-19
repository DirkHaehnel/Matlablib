function Z = FCSZfun(a,b,d,t1,t2,t3)

pp = sqrt(pi/2)^1.5*a*sqrt(b);
a = 4*d/a^2;
b = 4*d/b^2;
if nargin==3,
    Z = pp;
elseif nargin==4,
    Z = pp./2^1.5./(1+a*t1)./sqrt(1+b*t1);
elseif nargin==5,
    Z = pp./(3+4*a*(t1+t2)+4*a.^2*(t1.*t2))./sqrt(3+4*b*(t1+t2)+4*b.^2*(t1.*t2));
elseif nargin==6,
    Z = pp./8./(1+a*(3*t1+4*t2+3*t3)/2+4*a.^2*(t1.*t2+t2.*t3+t3.*t1)+2*a.^3*(t1.*t2.*t3))./...
        sqrt(1+b*(3*t1+4*t2+3*t3)/2+4*b.^2*(t1.*t2+t2.*t3+t3.*t1)+2*b.^3*(t1.*t2.*t3));
end


