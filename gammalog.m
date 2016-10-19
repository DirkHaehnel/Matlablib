function [f] = gammalog(z)
% GAMMALOG  Natural Log of the Gamma function valid in the entire complex plane.
%           This routine uses an excellent Lanczos series
%           for the complex Gamma function.
%
%usage: [f] = gammalog(z)
%             z may be complex and of any size.
%             Also  n! = prod(1:n) = exp(gammalog(n+1))
%        
%tested under version 5.3.1
%
%References: C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%            Y. Luke, "The Special ...", 1969 pp. 29-31
%            J. Spouge,  SIAM JNA 31, 1994. pp. 931
%            W. Press,  "Numerical Recipes"
%
%see also:   GAMMA GAMMALN GAMMAINC PSI
%see also:   mhelp GAMMA
%see also:   mhelp lnGAMMA

%Paul Godfrey
%pgodfrey@intersil.com
%7-28-00

twopi=pi+pi;
[row, col] = size(z);
z=z(:);
zz=z;

f = 0.*z; % reserve space in advance

p=find(real(z)<0);
if ~isempty(p)
   z(p)=-z(p);
end

%Lanczos approximation for the complex plane
c = [    0.9999999999995183;
       676.5203681218835;
     -1259.139216722289;
       771.3234287757674;
      -176.6150291498386;
        12.50734324009056;
        -0.1385710331296526;
         0.000009934937113930748;
         0.0000001659470187408462];

s=0;
t=z+7;
%Num Recipes used t=z+5
%for a lower order approximation

for k=9:-1:2
    s=s+c(k)./t;
    t=t-1;
end
s=s+c(1);
s=log(s) + log(pi+pi)/2 - (z+6.5) + (z-0.5).*log(z+6.5);

LogofGamma = s;
f = LogofGamma;

if ~isempty(p)
   f(p)=log(-pi)-log(zz(p))-f(p)-log(sin(pi*zz(p)));
end

 
p=find(round(zz)==zz & imag(zz)==0 & real(zz)<=0);
if ~isempty(p)
   f(p)=Inf;
end

p=find(real(zz)>=0 & imag(zz)==0);
if ~isempty(p)
pp=0;
for zint=1:172
    p=find(zz==zint);
    if ~isempty(p)
       f(p)=pp;
    end
    pp=pp+log(zint);
end
end

f=reshape(f,row,col);

return

%A demo of this routine is:
clc
clear all
close all
x=-4:1/16:4.5;
y=-4:1/16:4;
[X,Y]=meshgrid(x,y);
z=X+i*Y;
f=gammalog(z);
 
mesh(x,y,real(f),imag(f));
view([-30 15]);
rotate3d;
return
