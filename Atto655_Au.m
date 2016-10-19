[res0, head0] = tttrpro1('C:\Daten\Gregor\041104\atto655_01.t3r', [20 10]);
[res1, head1] = tttrpro1('C:\Daten\Gregor\041104\atto655_Au_5µw_01.t3r', [20 10]);
[res2, head2] = tttrpro1('C:\Daten\Gregor\041104\atto655_Au_2µw_01.t3r', [20 10]);

save atto655_Au

load atto655_Au

t0 = [1:45 50:70];
t1 = 1:70;
t2 = 1:12;

ind = res0.tau>116 & res0.tau<135;
tau = res0.tau(ind);
tst = sum(res0.tcspc1(ind,t0)')';
p01 = simplex('expfun',[3 1],[0 0],[],[],[],tau,tst-min(tst)); 
[err,c01] = expfun(p01,tau,tst-min(tst));
tst = sum(res0.tcspc2(ind,t0)')';
p02 = simplex('expfun',[3 1],[0 0],[],[],[],tau,tst-min(tst)); 
[err,c02] = expfun(p02,tau,tst-min(tst));

ind = res1.tau>116 & res1.tau<135;
tau = res1.tau(ind);
tst = sum(res1.tcspc1(ind,t1)')';
p11 = simplex('expfun',[3 1],[0 0],[],[],[],tau,tst-min(tst)); 
[err,c11] = expfun(p11,tau,tst-min(tst));
tst = sum(res1.tcspc2(ind,t1)')';
p12 = simplex('expfun',[3 1],[0 0],[],[],[],tau,tst-min(tst)); 
[err,c12] = expfun(p12,tau,tst-min(tst));

ind = res2.tau>116 & res2.tau<135;
tau = res2.tau(ind);
tst = sum(res2.tcspc1(ind,t2)')';
p21 = simplex('expfun',[3 1],[0 0],[],[],[],tau,tst-min(tst)); 
[err,c21] = expfun(p21,tau,tst-min(tst));
tst = sum(res2.tcspc2(ind,t2)')';
p22 = simplex('expfun',[3 1],[0 0],[],[],[],tau,tst-min(tst)); 
[err,c22] = expfun(p22,tau,tst-min(tst));

save atto655_Au t0 t1 t2 p01 p02 p11 p12 p21 p22 c01 c02 c11 c12 c21 c22 -append