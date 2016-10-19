name = '\\Mt200\MT200 data\MTusers\Mira\Area_16.pt3';

close all hidden

[intens, flim, tcspc, head] = PIScanRead(name); 
tim = FLIM2Mean(flim,tcspc); 
mim(mConv2(tim',Disk(1))*head.Resolution,log(1+intens'),[2.5 3.5],'v')
ylabel('time (ns)')

figure
mim(intens','v')

figure
t = head.Resolution*(0:2^12-1);
tau = Simplex('ExpFun',2,0,[],[],[],t(t>7.5 & t<20),tcspc(t>7.5 & t<20),1);
ax = axis;
text(ax(1)+0.8*diff(ax(1:2)),ax(3)+0.8*diff(ax(3:4)),[mnum2str(tau,1,1) ' ns'])
xlabel('time (ns)')
ylabel('intensity (a.u.)')

