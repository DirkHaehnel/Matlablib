load D:\Doc\FCS\Antibunching\Ingo4_2.mat

auto = res.auto(3:end,:);
bunch = anti.auto;
autotime = res.autotime(3:end);
bunchtime = 1e9*anti.autotime;

ind = bunchtime<=25+37.5; bunch = bunch(ind,:); bunchtime = bunchtime(ind);

close all

ind = autotime<=1e-5;
p = Simplex('ExpFun',1e-5,0,[],[],[],autotime(ind),sum(auto(ind,:),2),1);
for k=1:2
    p = Simplex('ExpFun',p,0*p,[],[],[],autotime(ind),sum(auto(ind,:),2),1);
end
off = mean(sum(auto(end-10:end,:),2));
tst=simplex('Rigler',[1e-2 1e-3 1./p'],[0 0 1./p'],[inf inf 1./p'],[],[],autotime,sum(auto,2),[],1);
for k=1:2
    tst = simplex('Rigler',tst,0*tst,[],[],[],autotime,sum(auto,2),[],1);
end
[err,c]=Rigler(tst,autotime,sum(auto,2),[],1);
[err,zz,zz]=Rigler(tst,c/c(1),autotime);
semilogx(autotime,sum(auto,2)/c(1),'o',autotime,cumsum(zz.*(ones(size(zz,1),1)*c')/c(1),2))
ax = axis;
axis([autotime(1) autotime(end) ax(3:4)])
xlabel('lag time [s]')
ylabel('autocorrelation')

pp = simplex('AntiBunch',[25 3 2],[],[],[],[],bunchtime,bunch,25,[],1);
for k=1:2
    pp = simplex('AntiBunch',pp,[],[],[],[],bunchtime,bunch,25,[],1);
end
[err,cc,c0,z,z,z] = AntiBunch(pp,bunchtime,bunch,25);
plot([-bunchtime(end:-1:2); bunchtime(1:end)],[bunch(end:-1:2,2); bunch(1:end,1)],'o',[-bunchtime(end:-1:1); bunchtime(2:end)],z(:,1))
ax = axis;
axis([pp(1)-37.5 pp(1)+37.5 ax(3:4)])
xlabel('time [ns]');
ylabel('autocorrelation [(cnts/s)^2]')

plot([-bunchtime(end:-1:2); bunchtime(1:end)],[bunch(end:-1:2,4); bunch(1:end,3)],'o',[-bunchtime(end:-1:1); bunchtime(2:end)],z(:,2))
axis([pp(1)-37.5 pp(1)+37.5 ax(3:4)])
xlabel('lag time \pm \Deltat [ns]');
ylabel('autocorrelation [(cnts/s)^2]')

lam = c(2)/sum(c(2:end))
h = lam.*(cc(2,:)-c0(2,:)*60/57)
delta = (sum(c(1:2))/sum(c)-1).*cc(2,:) + cc(3,:)
nemit = h/delta

