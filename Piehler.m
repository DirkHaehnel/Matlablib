% program Piehler

bin = 0:10;

res = pt3pro('D:\Daten\Gregor\Piehler_111705\R1N23_OG+TrisNTA_Atto_1.pt3',...
    'D:\Daten\Gregor\Piehler_111705\R1N23_OG+TrisNTA_Atto_2.pt3','bin');
h1 = mhist(res.bin1+res.bin2,bin);

res = pt3pro('D:\Daten\Gregor\Piehler_111705\R1N23_OG+Lig+TrisNTA_Atto_1.pt3',...
    'D:\Daten\Gregor\Piehler_111705\R1N23_OG+Lig+TrisNTA_Atto_2.pt3','bin');
h2 = mhist(res.bin1+res.bin2,bin);

res = pt3pro('D:\Daten\Gregor\Piehler_111705\R1N23_OG2_1.pt3',...
    'D:\Daten\Gregor\Piehler_111705\R1N23_OG2_2.pt3','bin');
h3 = mhist(res.bin1+res.bin2,bin);

clear res

res1 = pt3pro('D:\Daten\Gregor\Piehler_111705\R1N23_OG+TrisNTA_Atto_1.pt3',...
    'D:\Daten\Gregor\Piehler_111705\R1N23_OG+TrisNTA_Atto_2.pt3','xfcs');
res2 = pt3pro('D:\Daten\Gregor\Piehler_111705\R1N23_OG+Lig+TrisNTA_Atto_1.pt3',...
    'D:\Daten\Gregor\Piehler_111705\R1N23_OG+Lig+TrisNTA_Atto_2.pt3','xfcs');
res3 = pt3pro('D:\Daten\Gregor\Piehler_111705\R1N23_OG2_1.pt3',...
    'D:\Daten\Gregor\Piehler_111705\R1N23_OG2_2.pt3','xfcs');

save R1N23 h1 h2 h3 bin res1 res2 res3 


ind=res1.autotime>1e-5;
p1 = simplex('rigler',[1e-2 1e-3],[],[],[],[],res1.autotime(ind),sum(res1.auto(ind,:),2),[],1);
p2 = simplex('rigler',p1,[],[],[],[],res1.autotime(ind),sum(res2.auto(ind,:),2),[],1);
p3 = simplex('rigler',p1,[],[],[],[],res1.autotime(ind),sum(res3.auto(ind,:),2),[],1);

[err,b1] = rigler(p1,res1.autotime(ind),sum(res1.auto(ind,:),2),[],1);
[err,b2] = rigler(p2,res1.autotime(ind),sum(res2.auto(ind,:),2),[],1);
[err,b3] = rigler(p3,res1.autotime(ind),sum(res3.auto(ind,:),2),[],1);

ind=res1.tau>4;
a1 = simplex('expfun',[1 2],[],[],[],[],res1.tau(ind),sum(res1.tcspc(ind,:),2),1);
a2 = simplex('expfun',[1 2],[],[],[],[],res2.tau(ind),sum(res2.tcspc(ind,:),2),1);
a3 = simplex('expfun',[1 2],[],[],[],[],res3.tau(ind),sum(res3.tcspc(ind,:),2),1);

[err,c1] = expfun(a1,res1.tau(ind),sum(res1.tcspc(ind,:),2),1);
[err,c2] = expfun(a2,res2.tau(ind),sum(res2.tcspc(ind,:),2),1);
[err,c3] = expfun(a3,res3.tau(ind),sum(res3.tcspc(ind,:),2),1);

f1 = b1(2)/b1(1)*(sum(c1(2:3))/sum(c1))^2;
f2 = b2(2)/b2(1)*(sum(c2(2:3))/sum(c2))^2;
f3 = b3(2)/b3(1)*(sum(c3(2:3))/sum(c3))^2;

plot(bin/f1,h1/sum(h1),bin/f2,h2/sum(h2),bin/f3,h3/sum(h3))
xlabel('rel. brightness');
ylabel('probability density')
legend({'R1N23-OG+TrisNTA-Atto', 'R1N23-OG+Lig+TrisNTA-Atto', 'R1N23-OG'})


semilogy(res1.tau,res1.tcspc(:,1)/max(res1.tcspc(:,1)),res2.tau,res2.tcspc(:,1)/max(res2.tcspc(:,1)),res3.tau,res3.tcspc(:,1)/max(res3.tcspc(:,1)))
xlabel('time [ns]');
ylabel('normalized fluorecence intensity')
legend({'R1N23-OG+TrisNTA-Atto', 'R1N23-OG+Lig+TrisNTA-Atto', 'R1N23-OG'})
