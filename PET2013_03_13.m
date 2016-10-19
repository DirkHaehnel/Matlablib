[res, head] = SingleFocusCW2FCS('m:\MTusers\Qui\130313_FG_PET_cw\130313_FG-PET_cw_150mMNaCl_Fr03.ht3');

semilogx(1e-9*res.autotime1,sum(res.auto1(:,1,2,:)+res.auto1(:,2,1,:),4),res.autotime,sum(res.auto(:,1,2,:)+res.auto(:,2,1,:),4))
ax = axis;
axis([1e-10 1 0 ax(4)])
xlabel('time (s)'); 
ylabel('correlation')

ind = res.autotime>1e-8;

p = [1e-2 1e-3 1e6];
close; 
p = Simplex('Rigler',p,[],[],[],[],res.autotime(ind),sum(res.auto(ind,1,2,:)+res.auto(ind,2,1,:),4)-min(sum(res.auto(ind,1,2,:)+res.auto(ind,2,1,:),4)),[],1);
p = p([1 2 3 3]).*[1 1 1 0.5]';
for j=1:5
    p = Simplex('Rigler',p,[],[],[],[],res.autotime(ind),sum(res.auto(ind,1,2,:)+res.auto(ind,2,1,:),4)-min(sum(res.auto(ind,1,2,:)+res.auto(ind,2,1,:),4)),[],1);
end
[err, c, z] = Rigler(p,res.autotime(ind),sum(res.auto(ind,1,2,:)+res.auto(ind,2,1,:),4)-min(sum(res.auto(ind,1,2,:)+res.auto(ind,2,1,:),4)),[],1);
ax = axis;
axis([1e-8 1 0 ax(4)])
xlabel('time (s)');
ylabel ('correlation')
text(1e-4,0.9*ax(4),['\tau_1 = ' mnum2str(1e6/p(3),3,1) ' µs,      A_1 = ' mnum2str(c(3)/sum(c(2:4)),1,3)])
text(1e-4,0.8*ax(4),['\tau_2 = ' mnum2str(1e6/p(4),3,1) ' µs,      A_2 = ' mnum2str(c(4)/sum(c(2:4)),1,3)])



