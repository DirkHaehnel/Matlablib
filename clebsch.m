function cleb=clebsch(A,a,B,b,C,c)
%
%Clebsch-Gordan coefficients, see
%Brink & Satchler, P.140-->special case
%
%Brink & Satchler, P.34, eq.2.34-->general case
%
f1=0;
for k=max([0,B-C-a,A-C+b]):min([A+B-C,A-a,B+b])
  t1=prod(1:k);
  t2=prod(1:A+B-C-k);
  t3=prod(1:A-a-k);
  t4=prod(1:B+b-k); 
  t5=prod(1:C-B+a+k);
  t6=prod(1:C-A-b+k);

  f1=f1+(-1)^k*(t1*t2*t3*t4*t5*t6)^(-1);
end

d1=prod(1:-C+A+B);
d2=prod(1:C-A+B);
d3=prod(1:C+A-B);
d4=prod(1:1+C+A+B);

f2=d1*d2*d3/d4;

g1=prod(1:C-c);
g2=prod(1:C+c);
g3=prod(1:A-a);
g4=prod(1:A+a);
g5=prod(1:B-b);
g6=prod(1:B+b);

f3=(1+2*C)*g1*g2*g3*g4*g5*g6;

if c~=a+b || C<abs(A-B) || C>A+B || abs(a)>A || abs(b)>B || abs(c)>C
   cleb=0;
else
   cleb=f1*sqrt(f2*f3);
end;

if C==0 && c==0 
   if b==-a
      cleb=(-1)^(A-a)/sqrt(2*A+1);
   else
      cleb=0;
   end
end
