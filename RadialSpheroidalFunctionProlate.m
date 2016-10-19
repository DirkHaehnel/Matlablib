function [r1f,r1d,r2f,r2d] = RadialSpheroidalFunctionProlate(n,m,c,x,kf)

%       ==========================================================
%       Purpose: Compute oblate radial functions of the first
%       and second kinds, and their derivatives
%       Input :  m  --- Mode parameter,  m = 0,1,2,...
%       n  --- Mode parameter,  n = m,m+1,m+2,...
%       c  --- Spheroidal parameter
%       x  --- Argument (x ע 0)
%       cv --- Characteristic value
%       KF --- Function code
%       KF=1 for the first kind
%       KF=2 for the second kind
%       KF=3 for both the first and second kinds
%       Output:  R1F --- Radial function of the first kind
%       R1D --- Derivative of the radial function of
%       the first kind
%       R2F --- Radial function of the second kind
%       R2D --- Derivative of the radial function of
%       the second kind
%       Routines called:
%       (1) SDMN for computing expansion coefficients dk
%       (2) RMN1 for computing prolate or oblate radial
%       function of the first kind
%       (3) RMN2L for computing prolate or oblate radial
%       function of the second kind for a large argument
%       (4) RMN2SO for computing oblate radial functions of
%       the second kind for a small argument
%       ==========================================================

r1f=[]; r1d=[]; r2f=[]; r2d=[];
[cv,eg] = segv(m,n,c,-1);

df=[];id=[];
kd=1;
[m,n,c,cv,kd,df]=sdmn(m,n,c,cv,kd,df);
if (kf ~= 2);
    [m,n,c,x,df,kd,r1f,r1d]=rmn1(m,n,c,x,df,kd,r1f,r1d);
end;
if (kf > 1);
    [m,n,c,x,df,kd,r2f,r2d,id]=rmn2l(m,n,c,x,df,kd,r2f,r2d,id);
    if (id > -8);
        [m,n,c,x,cv,df,kd,r2f,r2d]=rmn2sp(m,n,c,x,cv,df,kd,r2f,r2d);
    end;
end;



function [m,n,c,cv,kd,df]=sdmn(m,n,c,cv,kd,df);

%     =====================================================
%     Purpose: Compute the expansion coefficients of the
%     prolate and oblate spheroidal functions, dk
%     Input :  m  --- Mode parameter
%     n  --- Mode parameter
%     c  --- Spheroidal parameter
%     cv --- Characteristic value
%     KD --- Function code
%     KD=1 for prolate; KD=-1 for oblate
%     Output:  DF(k) --- Expansion coefficients dk;
%     DF(1), DF(2), ... correspond to
%     d0, d2, ... for even n-m and d1,
%     d3, ... for odd n-m
%     =====================================================



sw=0.0;
fl=0.0;
nm=25+fix(0.5.*(n-m)+c);
if (c < 1.0d-10);
    for  i=1:nm;
        df(i)=0d0;
    end;
    df((n-m)./2+1)=1.0d0;
    return;
end;
cs=c.*c.*kd;
ip=1;
if (n-m == 2.*fix((n-m)./2)) ip=0; end;
for  i=1:nm+2;
    if (ip == 0) k=2.*(i-1); end;
    if (ip == 1) k=2.*i-1; end;
    dk0=m+k;
    dk1=m+k+1;
    dk2=2.*(m+k);
    d2k=2.*m+k;
    a(i)=(d2k+2.0).*(d2k+1.0)./((dk2+3.0).*(dk2+5.0)).*cs;
    d(i)=dk0.*dk1+(2.0.*dk0.*dk1-2.0.*m.*m-1.0)./((dk2-1.0) .*(dk2+3.0)).*cs;
    g(i)=k.*(k-1.0)./((dk2-3.0).*(dk2-1.0)).*cs;
end;
fs=1.0d0;
f1=0.0d0;
f0=1.0d-100;
kb=0;
df(nm+1)=0.0d0;
for  k=nm:-1:1;
    f=-((d(k+1)-cv).*f0+a(k+1).*f1)./g(k+1);
    if (abs(f) > abs(df(k+1)));
        df(k)=f;
        f1=f0;
        f0=f;
        if (abs(f) > 1.0d+100);
            for  k1=k:nm;
                df(k1)=df(k1).*1.0d-100;
            end;
            f1=f1.*1.0d-100;
            f0=f0.*1.0d-100;
        end;
    else;
        kb=k;
        fl=df(k+1);
        f1=1.0d-100;
        f2=-(d(1)-cv)./a(1).*f1;
        df(1)=f1;
        if (kb == 1);
            fs=f2;
        elseif (kb == 2);
            df(2)=f2;
            fs=-((d(2)-cv).*f2+g(2).*f1)./a(2);
        else;
            df(2)=f2;
            for  j=3:kb+1;
                f=-((d(j-1)-cv).*f2+g(j-1).*f1)./a(j-1);
                if (j <= kb) df(j)=f; end;
                if (abs(f) > 1.0d+100);
                    for  k1=1:j;
                        df(k1)=df(k1).*1.0d-100;
                    end;
                    f=f.*1.0d-100;
                    f2=f2.*1.0d-100;
                end;
                f1=f2;
                f2=f;
            end;
            fs=f;
        end;
        break;
    end;
end;
su1=0.0d0;
r1=1.0d0;
for  j=m+ip+1:2.*(m+ip);
    r1=r1.*j;
end;
su1=df(1).*r1;
for  k=2:kb;
    r1=-r1.*(k+m+ip-1.5d0)./(k-1.0d0);
    su1=su1+r1.*df(k);
end;
su2=0.0d0;
for  k=kb+1:nm;
    if (k ~= 1) r1=-r1.*(k+m+ip-1.5d0)./(k-1.0d0); end;
    su2=su2+r1.*df(k);
    if (abs(sw-su2) < abs(su2).*1.0d-14) break; end;
    sw=su2;
end;
r3=1.0d0;
for  j=1:(m+n+ip)./2;
    r3=r3.*(j+0.5d0.*(n+m+ip));
end;
r4=1.0d0;
for  j=1:(n-m-ip)./2;
    r4=-4.0d0.*r4.*j;
end;
s0=r3./(fl.*(su1./fs)+su2)./r4;
for  k=1:kb;
    df(k)=fl./fs.*s0.*df(k);
end;
for  k=kb+1:nm;
    df(k)=s0.*df(k);
end;



function [m,n,c,x,df,kd,r1f,r1d]=rmn1(m,n,c,x,df,kd,r1f,r1d);

%     =======================================================
%     Purpose: Compute prolate and oblate spheroidal radial
%     functions of the first kind for given m, n, c and x
%     Routines called:
%     (1) SCKB for computing expansion coefficients c2k
%     (2) SPHJ for computing the spherical Bessel
%     functions of the first kind
%     =======================================================



ck=[];sj=[];dj=[];
sw=0.0;
eps=1.0d-14;
ip=1;
nm1=fix((n-m)./2);
if (n-m == 2.*nm1) ip=0; end;
nm=25+nm1+fix(c);
reg=1.0d0;
if (m+nm > 80) reg=1.0d-200; end;
r0=reg;
for  j=1:2.*m+ip;
    r0=r0.*j;
end;
r=r0;
suc=r.*df(1);
for  k=2:nm;
    r=r.*(m+k-1.0).*(m+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
    suc=suc+r.*df(k);
    if (k > nm1&abs(suc-sw) < abs(suc).*eps) break; end;
    sw=suc;
end;
if (x == 0.0);
    [m,n,c,df,ck]=sckb(m,n,c,df,ck);
    sum=0.0d0;
    for  j=1:nm;
        sum=sum+ck(j);
        if (abs(sum-sw1) < abs(sum).*eps) break; end;
        sw1=sum;
    end;
    r1=1.0d0;
    for  j=1:(n+m+ip)./2;
        r1=r1.*(j+0.5d0.*(n+m+ip));
    end;
    r2=1.0d0;
    for  j=1:m;
        r2=2.0d0.*c.*r2.*j;
    end;
    r3=1.0d0;
    for  j=1:(n-m-ip)./2;
        r3=r3.*j;
    end;
    sa0=(2.0.*(m+ip)+1.0).*r1./(2.0.^n.*c.^ip.*r2.*r3);
    if (ip == 0);
        r1f=sum./(sa0.*suc).*df(1).*reg;
        r1d=0.0d0;
    elseif (ip == 1);
        r1f=0.0d0;
        r1d=sum./(sa0.*suc).*df(1).*reg;
    end;
    return;
end;
cx=c.*x;
nm2=2.*nm+m;
[nm2,cx,nm2,sj,dj]=sphj(nm2,cx,nm2,sj,dj);
a0=(1.0d0-kd./(x.*x)).^(0.5d0.*m)./suc;
r1f=0.0d0;
for  k=1:nm;
    l=2.*k+m-n-2+ip;
    if (l == 4.*fix(l./4)) lg=1; end;
    if (l ~= 4.*fix(l./4)) lg=-1; end;
    if (k == 1);
        r=r0;
    else;
        r=r.*(m+k-1.0).*(m+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
    end;
    np=m+2.*k-2+ip;
    r1f=r1f+lg.*r.*df(k).*sj(np+1);
    if (k > nm1&abs(r1f-sw) < abs(r1f).*eps) break; end;
    sw=r1f;
end;
r1f=r1f.*a0;
b0=kd.*m./x.^3.0d0./(1.0-kd./(x.*x)).*r1f;
sud=0.0d0;
for  k=1:nm;
    l=2.*k+m-n-2+ip;
    if (l == 4.*fix(l./4)) lg=1; end;
    if (l ~= 4.*fix(l./4)) lg=-1; end;
    if (k == 1);
        r=r0;
    else;
        r=r.*(m+k-1.0).*(m+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
    end;
    np=m+2.*k-2+ip;
    sud=sud+lg.*r.*df(k).*dj(np+1);
    if (k > nm1&abs(sud-sw) < abs(sud).*eps) break; end;
    sw=sud;
end;
r1d=b0+a0.*c.*sud;



function [m,n,c,df,ck]=sckb(m,n,c,df,ck);

%     ======================================================
%     Purpose: Compute the expansion coefficients of the
%     prolate and oblate spheroidal functions, c2k
%     Input :  m  --- Mode parameter
%     n  --- Mode parameter
%     c  --- Spheroidal parameter
%     DF(k) --- Expansion coefficients dk
%     Output:  CK(k) --- Expansion coefficients ck;
%     CK(1), CK(2), ... correspond to
%     c0, c2, ...
%     ======================================================




if (c <= 1.0d-10) c=1.0d-10; end;
nm=25+fix(0.5.*(n-m)+c);
ip=1;
if (n-m == 2.*fix((n-m)./2)) ip=0; end;
reg=1.0d0;
if (m+nm > 80) reg=1.0d-200; end;
fac=-0.5d0.^m;
for  k=0:nm-1;
    fac=-fac;
    i1=2.*k+ip+1;
    r=reg;
    for  i=i1:i1+2.*m-1;
        r=r.*i;
    end;
    i2=k+m+ip;
    for  i=i2:i2+k-1;
        r=r.*(i+0.5d0);
    end;
    sum=r.*df(k+1);
    for  i=k+1:nm;
        d1=2.0d0.*i+ip;
        d2=2.0d0.*m+d1;
        d3=i+m+ip-0.5d0;
        r=r.*d2.*(d2-1.0d0).*i.*(d3+k)./(d1.*(d1-1.0d0).*(i-k).*d3);
        sum=sum+r.*df(i+1);
        if (abs(sw-sum) < abs(sum).*1.0d-14) break; end;
        sw=sum;
    end;
    r1=reg;
    for  i=2:m+k;
        r1=r1.*i;
    end;
    ck(k+1)=fac.*sum./r1;
end;



function [n,x,nm,sj,dj]=sphj(n,x,nm,sj,dj);

%     =======================================================
%     Purpose: Compute spherical Bessel functions jn(x) and
%     their derivatives
%     Input :  x --- Argument of jn(x)
%     n --- Order of jn(x)  ( n = 0,1,תתת )
%     Output:  SJ(n) --- jn(x)
%     DJ(n) --- jn'(x)
%     NM --- Highest order computed
%     Routines called:
%     MSTA1 and MSTA2 for computing the starting
%     point for backward recurrence
%     =======================================================




nm=n;
if (abs(x) == 1.0d-100);
    for  k=0:n;
        sj(k+1)=0.0d0;
        dj(k+1)=0.0d0;
    end;
    sj(0+1)=1.0d0;
    dj(1+1)=.3333333333333333d0;
    return;
end;
sj(0+1)=sin(x)./x;
sj(1+1)=(sj(0+1)-cos(x))./x;
if (n >= 2);
    sa=sj(0+1);
    sb=sj(1+1);
    m=msta1(x,200);
    if (m < n);
        nm=m;
    else;
        m=msta2(x,n,15);
    end;
    f0=0.0d0;
    f1=1.0d0-100;
    for  k=m:-1:0;
        f=(2.0d0.*k+3.0d0).*f1./x-f0;
        if (k <= nm) sj(k+1)=f; end;
        f0=f1;
        f1=f;
    end;
    if (abs(sa) > abs(sb)) cs=sa./f; end;
    if (abs(sa) <= abs(sb)) cs=sb./f0; end;
    for  k=0:nm;
        sj(k+1)=cs.*sj(k+1);
    end;
end;
dj(0+1)=(cos(x)-sin(x)./x)./x;
for  k=1:nm;
    dj(k+1)=sj(k-1+1)-(k+1.0d0).*sj(k+1)./x;
end;



function [msta1]=msta1(x,mp);

%     ===================================================
%     Purpose: Determine the starting point for backward
%     recurrence such that the magnitude of
%     Jn(x) at that point is about 10^(-MP)
%     Input :  x     --- Argument of Jn(x)
%     MP    --- Value of magnitude
%     Output:  MSTA1 --- Starting point
%     ===================================================



a0=abs(x);
n0=fix(1.1.*a0)+1;
f0=envj(n0,a0)-mp;
n1=n0+5;
f1=envj(n1,a0)-mp;
for  it=1:20;
    nn=n1-(n1-n0)./(1.0d0-f0./f1);
    f=envj(nn,a0)-mp;
    if(abs(nn-n1) < 1) break; end;
    n0=n1;
    f0=f1;
    n1=nn;
    f1=f;
end;
msta1=fix(nn);



function [msta2]=msta2(x,n,mp);

%     ===================================================
%     Purpose: Determine the starting point for backward
%     recurrence such that all Jn(x) has MP
%     significant digits
%     Input :  x  --- Argument of Jn(x)
%     n  --- Order of Jn(x)
%     MP --- Significant digit
%     Output:  MSTA2 --- Starting point
%     ===================================================



a0=abs(x);
hmp=0.5d0.*mp;
ejn=envj(n,a0);
if (ejn <= hmp);
    obj=mp;
    n0=fix(1.1.*a0);
else;
    obj=hmp+ejn;
    n0=n;
end;
f0=envj(n0,a0)-obj;
n1=n0+5;
f1=envj(n1,a0)-obj;
for  it=1:20;
    nn=n1-(n1-n0)./(1.0d0-f0./f1);
    f=envj(nn,a0)-obj;
    if (abs(nn-n1) < 1) break; end;
    n0=n1;
    f0=f1;
    n1=nn;
    f1=f;
end;
msta2=fix(nn+10);


function [envj]=envj(n,x);



envj=0.5d0.*log10(6.28d0.*n)-n.*log10(1.36d0.*x./n);



function [m,n,c,x,df,kd,r2f,r2d,id]=rmn2l(m,n,c,x,df,kd,r2f,r2d,id);

%     ========================================================
%     Purpose: Compute prolate and oblate spheroidal radial
%     functions of the second kind for given m, n,
%     c and a large cx
%     Routine called:
%     SPHY for computing the spherical Bessel
%     functions of the second kind
%     ========================================================



sy=[];dy=[];
sw=0.0;
eps=1.0d-14;
ip=1;
nm1=fix((n-m)./2);
if (n-m == 2.*nm1) ip=0; end;
nm=25+nm1+fix(c);
reg=1.0d0;
if (m+nm > 80) reg=1.0d-200; end;
nm2=2.*nm+m;
cx=c.*x;
[nm2,cx,nm2,sy,dy]=sphy(nm2,cx,nm2,sy,dy);
r0=reg;
for  j=1:2.*m+ip;
    r0=r0.*j;
end;
r=r0;
suc=r.*df(1);
for  k=2:nm;
    r=r.*(m+k-1.0).*(m+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
    suc=suc+r.*df(k);
    if (k > nm1&abs(suc-sw) < abs(suc).*eps) break; end;
    sw=suc;
end;
a0=(1.0d0-kd./(x.*x)).^(0.5d0.*m)./suc;
r2f=0.0;
for  k=1:nm;
    l=2.*k+m-n-2+ip;
    if (l == 4.*fix(l./4)) lg=1; end;
    if (l ~= 4.*fix(l./4)) lg=-1; end;
    if (k == 1);
        r=r0;
    else;
        r=r.*(m+k-1.0).*(m+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
    end;
    np=m+2.*k-2+ip;
    r2f=r2f+lg.*r.*(df(k).*sy(np+1));
    eps1=abs(r2f-sw);
    if (k > nm1&eps1 < abs(r2f).*eps) break; end;
    sw=r2f;
end;
id1=fix(log10(eps1./abs(r2f)+eps));
r2f=r2f.*a0;
if (np >= nm2);
    id=10;
    return;
end;
b0=kd.*m./x.^3.0d0./(1.0-kd./(x.*x)).*r2f;
sud=0.0d0;
for  k=1:nm;
    l=2.*k+m-n-2+ip;
    if (l == 4.*fix(l./4)) lg=1; end;
    if (l ~= 4.*fix(l./4)) lg=-1; end;
    if (k == 1);
        r=r0;
    else;
        r=r.*(m+k-1.0).*(m+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
    end;
    np=m+2.*k-2+ip;
    sud=sud+lg.*r.*(df(k).*dy(np+1));
    eps2=abs(sud-sw);
    if (k > nm1&eps2 < abs(sud).*eps) break; end;
    sw=sud;
end;
r2d=b0+a0.*c.*sud;
id2=fix(log10(eps2./abs(sud)+eps));
id=max([id1,id2]);



function [n,x,nm,sy,dy]=sphy(n,x,nm,sy,dy);

%     ======================================================
%     Purpose: Compute spherical Bessel functions yn(x) and
%     their derivatives
%     Input :  x --- Argument of yn(x) ( x ע 0 )
%     n --- Order of yn(x) ( n = 0,1,תתת )
%     Output:  SY(n) --- yn(x)
%     DY(n) --- yn'(x)
%     NM --- Highest order computed
%     ======================================================




nm=n;
if (x < 1.0d-60);
    for  k=0:n;
        sy(k+1)=-1.0d+300;
        dy(k+1)=1.0d+300;
    end;
    return;
end;
sy(0+1)=-cos(x)./x;
sy(1+1)=(sy(0+1)-sin(x))./x;
f0=sy(0+1);
f1=sy(1+1);
for  k=2:n;
    f=(2.0d0.*k-1.0d0).*f1./x-f0;
    sy(k+1)=f;
    if (abs(f) >= 1.0d+300) break; end;
    f0=f1;
    f1=f;
end;
nm=k-1;
dy(0+1)=(sin(x)+cos(x)./x)./x;
for  k=1:nm;
    dy(k+1)=sy(k-1+1)-(k+1.0d0).*sy(k+1)./x;
end;



function [m,n,c,x,cv,df,kd,r2f,r2d]=rmn2sp(m,n,c,x,cv,df,kd,r2f,r2d);

%     ======================================================
%     Purpose: Compute prolate spheroidal radial function
%     of the second kind with a small argument
%     Routines called:
%     (1) LPMNS for computing the associated Legendre
%     functions of the first kind
%     (2) LQMNS for computing the associated Legendre
%     functions of the second kind
%     (3) KMN for computing expansion coefficients
%     and joining factors
%     ======================================================



dn=[];ck1=[];ck2=[];pm=[];pd=[];qm=[];qd=[];
if (abs(df(1)) < 1.0d-280);
    r2f=1.0d+300;
    r2d=1.0d+300;
    return;
end;
eps=1.0d-14;
ip=1;
nm1=fix((n-m)./2);
if (n-m == 2.*nm1) ip=0; end;
nm=25+nm1+fix(c);
nm2=2.*nm+m;
[m,n,c,cv,kd,df,dn,ck1,ck2]=kmn(m,n,c,cv,kd,df,dn,ck1,ck2);
[m,nm2,x,pm,pd]=lpmns(m,nm2,x,pm,pd);
[m,nm2,x,qm,qd]=lqmns(m,nm2,x,qm,qd);
su0=0.0d0;
for  k=1:nm;
    j=2.*k-2+m+ip;
    su0=su0+df(k).*qm(j+1);
    if (k > nm1&abs(su0-sw) < abs(su0).*eps) break; end;
    sw=su0;
end;
sd0=0.0d0;
for  k=1:nm;
    j=2.*k-2+m+ip;
    sd0=sd0+df(k).*qd(j+1);
    if (k > nm1&abs(sd0-sw) < abs(sd0).*eps) break; end;
    sw=sd0;
end;
su1=0.0d0;
sd1=0.0d0;
for  k=1:m;
    j=m-2.*k+ip;
    if (j < 0) j=-j-1; end;
    su1=su1+dn(k).*qm(j+1);
    sd1=sd1+dn(k).*qd(j+1);
end;
ga=((x-1.0d0)./(x+1.0d0)).^(0.5d0.*m);
for  k=1:m;
    j=m-2.*k+ip;
    if (~(j >= 0));
        if (j < 0) j=-j-1; end;
        r1=1.0d0;
        for  j1=1:j;
            r1=(m+j1).*r1;
        end;
        r2=1.0d0;
        for  j2=1:m-j-2;
            r2=j2.*r2;
        end;
        r3=1.0d0;
        sf=1.0d0;
        for  l1=1:j;
            r3=0.5d0.*r3.*(-j+l1-1.0).*(j+l1)./((m+l1).*l1).*(1.0-x);
            sf=sf+r3;
        end;
        if (m-j >= 2) gb=(m-j-1.0d0).*r2; end;
        if (m-j <= 1) gb=1.0d0; end;
        spl=r1.*ga.*gb.*sf;
        su1=su1+(-1).^(j+m).*dn(k).*spl;
        spd1=m./(x.*x-1.0d0).*spl;
        gc=0.5d0.*j.*(j+1.0)./(m+1.0);
        sd=1.0d0;
        r4=1.0d0;
        for  l1=1:j-1;
            r4=0.5d0.*r4.*(-j+l1).*(j+l1+1.0)./((m+l1+1.0).*l1) .*(1.0-x);
            sd=sd+r4;
        end;
        spd2=r1.*ga.*gb.*gc.*sd;
        sd1=sd1+(-1).^(j+m).*dn(k).*(spd1+spd2);
    end;
end;
su2=0.0d0;
ki=(2.*m+1+ip)./2;
nm3=nm+ki;
for  k=ki:nm3;
    j=2.*k-1-m-ip;
    su2=su2+dn(k).*pm(j+1);
    if (j > m&abs(su2-sw) < abs(su2).*eps) break; end;
    sw=su2;
end;
sd2=0.0d0;
for  k=ki:nm3;
    j=2.*k-1-m-ip;
    sd2=sd2+dn(k).*pd(j+1);
    if (j > m&abs(sd2-sw) < abs(sd2).*eps) break; end;
    sw=sd2;
end;
sum=su0+su1+su2;
sdm=sd0+sd1+sd2;
r2f=sum./ck2;
r2d=sdm./ck2;



function [m,n,x,pm,pd]=lpmns(m,n,x,pm,pd);

%     ========================================================
%     Purpose: Compute associated Legendre functions Pmn(x)
%     and Pmn'(x) for a given order
%     Input :  x --- Argument of Pmn(x)
%     m --- Order of Pmn(x),  m = 0,1,2,...,n
%     n --- Degree of Pmn(x), n = 0,1,2,...,N
%     Output:  PM(n) --- Pmn(x)
%     PD(n) --- Pmn'(x)
%     ========================================================




for  k=0:n;
    pm(k+1)=0.0d0;
    pd(k+1)=0.0d0;
end;
if (abs(x) == 1.0d0);
    for  k=0:n;
        if (m == 0);
            pm(k+1)=1.0d0;
            pd(k+1)=0.5d0.*k.*(k+1.0);
            if (x < 0.0);
                pm(k+1)=(-1).^k.*pm(k+1);
                pd(k+1)=(-1).^(k+1).*pd(k+1);
            end;
        elseif (m == 1);
            pd(k+1)=1.0d+300;
        elseif (m == 2);
            pd(k+1)=-0.25d0.*(k+2.0).*(k+1.0).*k.*(k-1.0);
            if (x < 0.0) pd(k+1)=(-1).^(k+1).*pd(k+1); end;
        end;
    end;
    return;
end;
x0=abs(1.0d0-x.*x);
pm0=1.0d0;
pmk=pm0;
for  k=1:m;
    pmk=(2.0d0.*k-1.0d0).*sqrt(x0).*pm0;
    pm0=pmk;
end;
pm1=(2.0d0.*m+1.0d0).*x.*pm0;
pm(m+1)=pmk;
pm(m+1+1)=pm1;
for  k=m+2:n;
    pm2=((2.0d0.*k-1.0d0).*x.*pm1-(k+m-1.0d0).*pmk)./(k-m);
    pm(k+1)=pm2;
    pmk=pm1;
    pm1=pm2;
end;
pd(0+1)=((1.0d0-m).*pm(1+1)-x.*pm(0+1))./(x.*x-1.0);
for  k=1:n;
    pd(k+1)=(k.*x.*pm(k+1)-(k+m).*pm(k-1+1))./(x.*x-1.0d0);
end;



function [m,n,x,qm,qd]=lqmns(m,n,x,qm,qd);

%     ========================================================
%     Purpose: Compute associated Legendre functions Qmn(x)
%     and Qmn'(x) for a given order
%     Input :  x --- Argument of Qmn(x)
%     m --- Order of Qmn(x),  m = 0,1,2,...
%     n --- Degree of Qmn(x), n = 0,1,2,...
%     Output:  QM(n) --- Qmn(x)
%     QD(n) --- Qmn'(x)
%     ========================================================




for  k=0:n;
    qm(k+1)=0.0d0;
    qd(k+1)=0.0d0;
end;
if (abs(x) == 1.0d0);
    for  k=0:n;
        qm(k+1)=1.0d+300;
        qd(k+1)=1.0d+300;
    end;
    return;
end;
ls=1;
if (abs(x) > 1.0d0) ls=-1; end;
xq=sqrt(ls.*(1.0d0-x.*x));
q0=0.5d0.*log(abs((x+1.0)./(x-1.0)));
q00=q0;
q10=-1.0d0./xq;
q01=x.*q0-1.0d0;
q11=-ls.*xq.*(q0+x./(1.0d0-x.*x));
qf0=q00;
qf1=q10;
for  k=2:m;
    qm0=-2.0d0.*(k-1.0)./xq.*x.*qf1-ls.*(k-1.0).*(2.0-k).*qf0;
    qf0=qf1;
    qf1=qm0;
end;
if (m == 0) qm0=q00; end;
if (m == 1) qm0=q10; end;
qm(0+1)=qm0;
if (abs(x) < 1.0001d0);
    if (m == 0&n > 0);
        qf0=q00;
        qf1=q01;
        for  k=2:n;
            qf2=((2.0.*k-1.0d0).*x.*qf1-(k-1.0).*qf0)./k;
            qm(k+1)=qf2;
            qf0=qf1;
            qf1=qf2;
        end;
    end;
    qg0=q01;
    qg1=q11;
    for  k=2:m;
        qm1=-2.0d0.*(k-1.0)./xq.*x.*qg1-ls.*k.*(3.0-k).*qg0;
        qg0=qg1;
        qg1=qm1;
    end;
    if (m == 0) qm1=q01; end;
    if (m == 1) qm1=q11; end;
    qm(1+1)=qm1;
    if (m == 1&n > 1);
        qh0=q10;
        qh1=q11;
        for  k=2:n;
            qh2=((2.0.*k-1.0d0).*x.*qh1-k.*qh0)./(k-1.0);
            qm(k+1)=qh2;
            qh0=qh1;
            qh1=qh2;
        end;
    elseif (m >= 2);
        qg0=q00;
        qg1=q01;
        qh0=q10;
        qh1=q11;
        for  l=2:n;
            q0l=((2.0d0.*l-1.0d0).*x.*qg1-(l-1.0d0).*qg0)./l;
            q1l=((2.0.*l-1.0d0).*x.*qh1-l.*qh0)./(l-1.0d0);
            qf0=q0l;
            qf1=q1l;
            for  k=2:m;
                qmk=-2.0d0.*(k-1.0)./xq.*x.*qf1-ls.*(k+l-1.0).* (l+2.0-k).*qf0;
                qf0=qf1;
                qf1=qmk;
            end;
            qm(l+1)=qmk;
            qg0=qg1;
            qg1=q0l;
            qh0=qh1;
            qh1=q1l;
        end;
    end;
else;
    if (abs(x) > 1.1);
        km=40+m+n;
    else;
        km=(40+m+n).*fix(-1.0-1.8.*log(x-1.0));
    end;
    qf2=0.0d0;
    qf1=1.0d0;
    for  k=km:-1:0;
        qf0=((2.0.*k+3.0d0).*x.*qf1-(k+2.0-m).*qf2)./(k+m+1.0);
        if (k <= n) qm(k+1)=qf0; end;
        qf2=qf1;
        qf1=qf0;
    end;
    for  k=0:n;
        qm(k+1)=qm(k+1).*qm0./qf0;
    end;
end;
if (abs(x) < 1.0d0);
    for  k=0:n;
        qm(k+1)=(-1).^m.*qm(k+1);
    end;
end;
qd(0+1)=((1.0d0-m).*qm(1+1)-x.*qm(0+1))./(x.*x-1.0);
for  k=1:n;
    qd(k+1)=(k.*x.*qm(k+1)-(k+m).*qm(k-1+1))./(x.*x-1.0);
end;



function [m,n,c,cv,kd,df,dn,ck1,ck2]=kmn(m,n,c,cv,kd,df,dn,ck1,ck2);

%     ===================================================
%     Purpose: Compute the expansion coefficients of the
%     prolate and oblate spheroidal functions
%     and joining factors
%     ===================================================



nm=25+fix(0.5.*(n-m)+c);
nn=nm+m;
cs=c.*c.*kd;
ip=1;
if (n-m == 2.*fix((n-m)./2)) ip=0; end;
for  i=1:nn+3;
    if (ip == 0) k=-2.*(i-1); end;
    if (ip == 1) k=-(2.*i-3); end;
    gk0=2.0d0.*m+k;
    gk1=(m+k).*(m+k+1.0d0);
    gk2=2.0d0.*(m+k)-1.0d0;
    gk3=2.0d0.*(m+k)+3.0d0;
    u(i)=gk0.*(gk0-1.0d0).*cs./(gk2.*(gk2+2.0d0));
    v(i)=gk1-cv+(2.0d0.*(gk1-m.*m)-1.0d0).*cs./(gk2.*gk3);
    w(i)=(k+1.0d0).*(k+2.0d0).*cs./((gk2+2.0d0).*gk3);
end;
for  k=1:m;
    t=v(m+1);
    for  l=0:m-k-1;
        t=v(m-l)-w(m-l+1).*u(m-l)./t;
    end;
    rk(k)=-u(k)./t;
end;
r=1.0d0;
for  k=1:m;
    r=r.*rk(k);
    dn(k)=df(1).*r;
end;
tp(nn)=v(nn+1);
for  k=nn-1:-1:m+1;
    tp(k)=v(k+1)-w(k+2).*u(k+1)./tp(k+1);
    if (k > m+1) rk(k)=-u(k)./tp(k); end;
end;
if (m == 0) dnp=df(1); end;
if (m ~= 0) dnp=dn(m); end;
dn(m+1)=(-1).^ip.*dnp.*cs./((2.0.*m-1.0).*(2.0.*m+1.0-4.0.*ip) .*tp(m+1));
for  k=m+2:nn;
    dn(k)=rk(k).*dn(k-1);
end;
r1=1.0d0;
for  j=1:(n+m+ip)./2;
    r1=r1.*(j+0.5d0.*(n+m+ip));
end;
nm1=(n-m)./2;
r=1.0d0;
for  j=1:2.*m+ip;
    r=r.*j;
end;
su0=r.*df(1);
for  k=2:nm;
    r=r.*(m+k-1.0).*(m+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
    su0=su0+r.*df(k);
    if (k > nm1&abs((su0-sw)./su0) < 1.0d-14) break; end;
    sw=su0;
end;
if (~(kd == 1));
    r2=1.0d0;
    for  j=1:m;
        r2=2.0d0.*c.*r2.*j;
    end;
    r3=1.0d0;
    for  j=1:(n-m-ip)./2;
        r3=r3.*j;
    end;
    sa0=(2.0.*(m+ip)+1.0).*r1./(2.0.^n.*c.^ip.*r2.*r3.*df(1));
    ck1=sa0.*su0;
    if (kd == -1) return; end;
end;
r4=1.0d0;
for  j=1:(n-m-ip)./2;
    r4=4.0d0.*r4.*j;
end;
r5=1.0d0;
for  j=1:m;
    r5=r5.*(j+m)./c;
end;
g0=dn(m);
if (m == 0) g0=df(1); end;
sb0=(ip+1.0).*c.^(ip+1)./(2.0.*ip.*(m-2.0)+1.0)./(2.0.*m-1.0);
ck2=(-1).^ip.*sb0.*r4.*r5.*g0./r1.*su0;



function [m,n,c,cv,eg]=segv(m,n,c,kd,cv,eg);

%     =========================================================
%     Purpose: Compute the characteristic values of spheroidal
%     wave functions
%     Input :  m  --- Mode parameter
%     n  --- Mode parameter
%     c  --- Spheroidal parameter
%     KD --- Function code
%     KD=1 for Prolate; KD=-1 for Oblate
%     Output:  CV --- Characteristic value for given m, n and c
%     EG(L) --- Characteristic value for mode m and n'
%     ( L = n' - m + 1 )
%     =========================================================



if (c < 1.0d-10);
    for  i=1:n;
        eg(i)=(i+m).*(i+m-1.0d0);
    end;
else;
    icm=fix((n-m+2)./2);
    nm=10+fix(0.5.*(n-m)+c);
    cs=c.*c.*kd;
    for  l=0:1;
        for  i=1:nm;
            if (l == 0) k=2.*(i-1); end;
            if (l == 1) k=2.*i-1; end;
            dk0=m+k;
            dk1=m+k+1;
            dk2=2.*(m+k);
            d2k=2.*m+k;
            a(i)=(d2k+2.0).*(d2k+1.0)./((dk2+3.0).*(dk2+5.0)).*cs;
            d(i)=dk0.*dk1+(2.0.*dk0.*dk1-2.0.*m.*m-1.0)./((dk2-1.0) .*(dk2+3.0)).*cs;
            g(i)=k.*(k-1.0)./((dk2-3.0).*(dk2-1.0)).*cs;
        end;
        for  k=2:nm;
            e(k)=sqrt(a(k-1).*g(k));
            f(k)=e(k).*e(k);
        end;
        f(1)=0.0d0;
        e(1)=0.0d0;
        xa=d(nm)+abs(e(nm));
        xb=d(nm)-abs(e(nm));
        nm1=nm-1;
        for  i=1:nm1;
            t=abs(e(i))+abs(e(i+1));
            t1=d(i)+t;
            if (xa < t1) xa=t1; end;
            t1=d(i)-t;
            if (t1 < xb) xb=t1; end;
        end;
        for  i=1:icm;
            b(i)=xa;
            h(i)=xb;
        end;
        for  k=1:icm;
            for  k1=k:icm;
                if (b(k1) < b(k));
                    b(k)=b(k1);
                    break;
                end;
            end;
            if (k ~= 1&h(k) < h(k-1)) h(k)=h(k-1); end;
            while (1==1);
                x1=(b(k)+h(k))./2.0d0;
                cv0(k)=x1;
                if (abs((b(k)-h(k))./x1) < 1.0d-14) break; end;
                j=0;
                s=1.0d0;
                for  i=1:nm;
                    if (s == 0.0d0) s=s+1.0d-30; end;
                    t=f(i)./s;
                    s=d(i)-t-x1;
                    if (s < 0.0d0) j=j+1; end;
                end;
                if (j < k);
                    h(k)=x1;
                else;
                    b(k)=x1;
                    if (j >= icm);
                        b(icm)=x1;
                    else;
                        if (h(j+1) < x1) h(j+1)=x1; end;
                        if (x1 < b(j)) b(j)=x1; end;
                    end;
                end;
            end;
            cv0(k)=x1;
            if (l == 0) eg(2.*k-1)=cv0(k); end;
            if (l == 1) eg(2.*k)=cv0(k); end;
        end;
    end;
end;
cv=eg(n-m+1);

