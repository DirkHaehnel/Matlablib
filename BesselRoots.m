function [rj0, rj1, ry0, ry1]  =  BesselRoots(n,nt)

%                Yn(x), and their derivatives
%       Input :  n  --- Order of Bessel functions ( n ó 101 )
%                NT --- Number of zeros (roots)
%       Output:  RJ0(L) --- L-th zero of Jn(x),  L = 1,2,...,NT
%                RJ1(L) --- L-th zero of Jn'(x), L = 1,2,...,NT
%                RY0(L) --- L-th zero of Yn(x),  L = 1,2,...,NT
%                RY1(L) --- L-th zero of Yn'(x), L = 1,2,...,NT
%       Routine called: JYNDD for computing Jn(x), Yn(x), and
%                       their first and second derivatives

x = []; bjn = []; djn = []; fjn = []; byn = []; dyn = []; fyn = [];
if n<=20
    x = 2.82141+1.15859.*n;
else;
    x = n+1.85576.*n.^0.33333+1.03315./n.^0.33333;
end;
l = 0;
while 1
    x0 = x;
    [n,x,bjn,djn,fjn,byn,dyn,fyn] = jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn);
    x = x-bjn./djn;
    if ~(abs(x-x0)>1e-9) break; end;
    l = l+1;
    rj0(l) = x;
    x = x+3.1416+(0.0972+0.0679.*n-0.000354.*n.^2)./l;
    if (~(l < nt)) break; end;
end;
if n<=20
    x = 0.961587+1.07703.*n;
else;
    x = n+0.80861.*n.^(1/3)+0.07249./n.^(1/3);
end;
if n==0 x = 3.8317; end;
l = 0;
while 1
    x0 = x;
    [n,x,bjn,djn,fjn,byn,dyn,fyn] = jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn);
    x = x-djn./fjn;
    if (~(abs(x-x0) > 1.0d-9)) break; end;
    l = l+1;
    rj1(l) = x;
    x = x+pi+(0.4955+0.0915.*n-0.000435.*n.^2)./l;
    if (~(l < nt)) break; end;
end;
if n<=20
    x = 1.19477+1.08933.*n;
else;
    x = n+0.93158.*n.^(1/3)+0.26035./n.^(1/3);
end;
l = 0;
while 1
    x0 = x;
    [n,x,bjn,djn,fjn,byn,dyn,fyn] = jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn);
    x = x-byn./dyn;
    if (~(abs(x-x0) > 1.0d-9)) break; end;
    l = l+1;
    ry0(l) = x;
    x = x+3.1416+(0.312+0.0852.*n-0.000403.*n.^2)./l;
    if ~(l<nt) break; end;
end;
if n<=20
    x = 2.67257+1.16099.*n;
else;
    x = n+1.8211.*n.^(1/3)+0.94001./n.^(1/3);
end;
l = 0;
while 1
    x0 = x;
    [n,x,bjn,djn,fjn,byn,dyn,fyn] = jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn);
    x = x-dyn./fyn;
    if ~(abs(x-x0)>1e-9) break; end;
    l = l+1;
    ry1(l) = x;
    x = x+pi+(0.197+0.0643.*n-0.000286.*n.^2)./l;
    if ~(l<nt) break; end;
end;
return;



function [n,x,bjn,djn,fjn,byn,dyn,fyn] = jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn);

%                their first and second derivatives
%       Input:   x   ---  Argument of Jn(x) and Yn(x) ( x > 0 )
%                n   ---  Order of Jn(x) and Yn(x)
%       Output:  BJN ---  Jn(x)
%                DJN ---  Jn'(x)
%                FJN ---  Jn"(x)
%                BYN ---  Yn(x)
%                DYN ---  Yn'(x)
%                FYN ---  Yn"(x)

bj = zeros(102,1);by = zeros(102,1);

for  nt = 1:900;
    mt = fix(0.5.*log10(6.28.*nt)-nt.*log10(1.36.*abs(x)./nt));
    if mt>20 break; end;
end;
m = nt;
bs = 0;
f0 = 0;
f1 = 1e-35;
su = 0;
for  k = m:-1:0;
    f = 2.*(k+1).*f1./x-f0;
    if k<=n+1 bj(k+1) = f; end;
    if k==2.*fix(k./2)
        bs = bs+2.*f;
        if k~=0 su = su+(-1).^(k./2).*f./k; end;
    end
    f0 = f1;
    f1 = f;
end
for  k = 0:n+1;
    bj(k+1) = bj(k+1)./(bs-f);
end
bjn = bj(n+1);
ec = 0.5772156649015329;
e0 = 0.3183098861837907;
s1 = 2.*e0.*(log(x./2)+ec).*bj(1);
f0 = s1-8.*e0.*su./(bs-f);
f1 = (bj(2).*f0-2.*e0./x)./bj(1);
by(1) = f0;
by(2) = f1;
for  k = 2:n+1;
    f = 2.*(k-1).*f1./x-f0;
    by(k+1) = f;
    f0 = f1;
    f1 = f;
end;
byn = by(n+1);
djn = -bj(n+2)+n.*bj(n+1)./x;
dyn = -by(n+2)+n.*by(n+1)./x;
fjn = (n.*n./(x.*x)-1).*bjn-djn./x;
fyn = (n.*n./(x.*x)-1).*byn-dyn./x;
return;

