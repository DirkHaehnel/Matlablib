function f = erfz(z)
%ERFZ Complex Error function
%     f = erfz(z) is the error function
%     for the elements of z.
%     Z may be complex and of any size.
%
%Usage:   f = erfz(z)
%
%Ref: Abramowitz & Stegun section 7.1
%
%Tested under version 5.2.0
%
%See also erf, erfc, erfcx, erfinc, erfcore


%Paul Godfrey
%pjg@mlb.semi.harris.com


twopi=2*pi;
[row col]=size(z);
z=z(:);
f=zeros(size(z));
nn=25;

x=real(z);
y=imag(z);
k1=2*exp(-x.*x)/pi;
k2=exp(-i*2*x.*y);

s1=erf(x);
s2=k1./(4*x).*(1-k2);
p=find(isnan(s2));
s2(p)=i*y(p)/pi;
s5=0;
for n=1:nn;
    s3=exp(-n*n/4)./(n*n+4*x.*x);
    s4=2*x-k2.*(2*x.*cosh(n*y)-i*n*sinh(n*y));
    s5=s5+s3.*s4;
end
s6=k1.*s5;
f=s1+s2+s6;

f=reshape(f,row,col);

return

%a demo of this function is
x=-4:0.125:4;
y=x;
[X,Y]=meshgrid(x,y);
z=X+i*Y;
f=erfz(z);
af=abs(f);
%let's truncate for visibility
p=find(af>5);
af(p)=5;
mesh(af);
axis([0 80 0 80 0 5]);
view(-70, 40);
rotate3d;

