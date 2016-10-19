% BuschmannMaster

if 0
    path = 'D:\Joerg\Doc\Hell\Buschmann\';
    names = dir([path '*.pt3']);

    for j=1:length(names)
        [res(j), head(j)] = Buschmann([path names(j).name], 'xfcs');
    end

    for j=1:length(names) auto(:,j)=res(j).auto; end
    t=res(1).autotime;
    for j=1:length(names) tcspc(:,j)=res(j).tcspc(:,1); end
    tau=res(1).tau;
    for j=1:5 p(:,j)=simplex('ExpFun',[2],[],[],[],[],tau,tcspc,1); end
    
    eval(['save ' path 'data res head p tau t tcspc auto names'])
end

if 1
    path = 'D:\Joerg\Doc\Hell\Buschmann\';
    % eval(['load ' path 'data'])
    ind=tau>3;   
    k = 4;
    [err, c, zz, z] = ExpFun(p(k), tau(ind), tcspc(ind,k),1);
    z1 = ExpFun(p(k), c.*[0 1]', tau(ind));
    z2 = ExpFun(p(k)/2, c.*[0 1]', tau(ind));
    
    k = 2;
    [res, head] = Buschmann([path names(k).name], 'flcs', [tau(ind)' z1 z2]);
    
end