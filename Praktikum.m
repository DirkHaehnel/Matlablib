cd d:/joerg/doc/gregor/diffstandards

[y, flag, x] = tttrread('cy5_900_150µW.t3r');

plot(y(1:1e3));
z = tttr2bin(y,1e4);

close
h=plot(z(1:1e3)); 
set(gcf,'renderer','zbuffer'); 
for j=1:1000
    set(h,'ydata',z(1+j*2e2:1e3+j*2e2)); 
    ax=axis; 
    axis([ax(1:2) 0 35]); 
    drawnow;  
end


mhist(x,min(x):max(x))
[x1, tc] = mhist(x(flag==0),min(x):max(x));
[x2, tc] = mhist(x(flag==1),min(x):max(x));
plot(tc,x1,tc,x2)


Ncasc = 17; Nsub = 10;
autotime = cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
auto=tttr2fcs(y, flag==1, flag==1, Ncasc, Nsub, 1);
auto1=tttr2fcs(y, flag==0, flag==1, Ncasc, Nsub, 1);
FCSScalePlot(autotime,[auto auto1])

load optsatdata
FCSScalePlot(xdata(15).autotime,[xdata(15).auto(:,1,2) xdata(16).auto(:,1,2) xdata(17).auto(:,1,2) xdata(18).auto(:,1,2)])

load data
tm1=(fluo.auto(10:end,1,2)-mean(fluo.auto(end-10:end,1,2)))/(fluo.auto(10,1,2)-mean(fluo.auto(end-10:end,1,2)));
t = fluo.autotime(10:end);
tt = t(1:end-30);
semilogx(t,tm1)

td1 = cy5.auto(10:end-30,1,2);
close; p=simplex('affineexpfit',[1e-2 1e5],[0 0],[],[],[],tt,td1,t,tm1)


