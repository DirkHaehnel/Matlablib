clear p1 p2 err1 err2
dv = 3:0.2:5; 
for j=1:length(dv) 
    drot=dv(j); 
    for k=1:10 
        p1(:,k,j)=Simplex('ExpFun',[2 drot],[0 drot],[inf drot],[],[],res.tau(ind),res.tcspc(ind,1)); 
        p2(:,k,j)=Simplex('ExpFun',[2 drot],[0 drot],[inf drot],[],[],res.tau(ind),res.tcspc(ind,2)); 
        err1(k,j)=ExpFun(p1(:,k,j),res.tau(ind),res.tcspc(ind,1)); 
        err2(k,j)=ExpFun(p2(:,k,j),res.tau(ind),res.tcspc(ind,2)); 
    end; 
end