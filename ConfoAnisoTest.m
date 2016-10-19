exc.NA = 0.001;
exc.n0 = 1;
exc.n = 1;
exc.n1 = 1;
exc.d0 = [];
exc.d = 0;
exc.d1 = [];
exc.focpos = 0;
exc.maxm = 0;
exc.rho = [0 0.01]';
exc.z = 0;
exc.fxc = [1 1]';
exc.fxs = [0 0]';
exc.fyc  = [0 0]';
exc.fys = [0 0]';
exc.fzc = [0 0]';
exc.fzs = [0 0]';
exc.maxnum = 1;

[ux0, uy0, pp0, po0, op0, oo0] = GaussExc2RotoFCS(exc, 0.5, 60, 100, 0);

po0 = po0/pp0(1,1);
op0 = op0/pp0(1,1);
oo0 = oo0/pp0(1,1);
pp0 = pp0/pp0(1,1);

t = (0:1e3)';
tst = RotoDiffCWFit1(40,t,[], 0, inf, 1);
tst1 = RotoDiffCWFit1(40,t,[], 0, 0.6, 1);
dp = 1/6/Rad2RotoDiff(40,20);
dd = 0;
zz = zeros(length(t),4);
for L=1:4
    for M=0:L
        tmp = exp(-L*(L+1)*dp*t-M^2*dd*t);
        zz = zz + [pp0(L+1,M+1)*tmp, po0(L+1,M+1)*tmp, op0(L+1,M+1)*tmp oo0(L+1,M+1)*tmp];
    end
end
semilogx(t,zz,t,tst/tst(1,1)*zz(1,1),'o',t,tst1/tst1(1,1)*zz(1,1),'x')

