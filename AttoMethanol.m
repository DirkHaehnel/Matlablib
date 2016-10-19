% AttoMethanol

load D:\Daten\Dertinger\AttoMethanol090505\atto655_meth.mat res

autotime=res.autotime;
auto=[res.auto(:,1,2)+res.auto(:,2,1) res.auto(:,3,4)+res.auto(:,4,3) (res.auto(:,1,4)+res.auto(:,4,1)+res.auto(:,2,3)+res.auto(:,3,2))/2];
tt=autotime(autotime>=1e-4);

ind = autotime<tt(end); 
p = simplex('MultiFocusCrossFit',[1e-3 1e-2 0 5e-3],[0 0 0 0],[inf inf 0 inf],[],[],autotime(ind),auto(ind,1),auto(ind,2),auto(ind,3));
for j=1:5
    p = simplex('MultiFocusCrossFit',p,[0 0 0 0],[inf inf 0 inf],[],[],autotime(ind),auto(ind,1),auto(ind,2),auto(ind,3));
end

for j=2:length(tt) 
    ind=autotime<tt(length(tt)-j+1); 
    p(:,j)=simplex('MultiFocusCrossFit',p(:,j-1),[0 0 0 0],[inf inf 0 inf],[],[],autotime(ind),auto(ind,1),auto(ind,2),auto(ind,3)); 
end