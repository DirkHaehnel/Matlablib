% BSAAnisotropy

load 'c:\Joerg\Doc\Fcs\BSA\BSA_anysotropy_194.mat'

bin = res.bin';
tcspcdata = res.tcspc;
tau = []; tcspc = [];
for d=1:2
    [ind, num] = mCluster(tcspcdata(:,d)-min(tcspcdata(:,d))>0.01*(max(tcspcdata(:,d))-min(tcspcdata(:,d))));
    tcspc(1:num(end),2*d-1) = tcspcdata(ind==length(num), d);
    tcspc(1:num(end-1),2*d) = tcspcdata(ind==length(num)-1, d);
    tau(1:num(end),2*d-1) = bin(ind==length(num))';
    tau(1:num(end-1),2*d) = bin(ind==length(num)-1)';
    if max(tau(:,2*d-1))>max(tau(:,2*d))
        tmp = tau(:,2*d-1); tau(:,2*d-1) = tau(:,2*d); tau(:,2*d) = tmp;
        tmp = tcspc(:,2*d-1); tcspc(:,2*d-1) = tcspc(:,2*d); tcspc(:,2*d) = tmp;        
    end
end
tcspc(any(tau==0,2),:)=[];
tau(any(tau==0,2),:)=[];
for j=1:4
    tcspc(:,j) = tcspc(:,j)-min(tcspc(:,j));
end
semilogy(tau,tcspc); drawnow

ind = ceil(2/head.Resolution);
p = Simplex('ExpFun',[22 40],[0 0],[],[],[],head.Resolution*tau(ind:end,1),sqrt(tcspc(ind:end,1:2).*tcspc(ind:end,4:-1:3)),1);
[err,c,zz,z] = ExpFun(p,head.Resolution*tau(ind:end,1),sqrt(tcspc(ind:end,1:2).*tcspc(ind:end,4:-1:3)));
plot(head.Resolution*tau(ind:end,1),sqrt(tcspc(ind:end,1:2).*tcspc(ind:end,4:-1:3)),'o',head.Resolution*tau(ind:end,1),z)
xlabel('time (ns)');
ylabel('intensity (a.u.)');
legend({'parallel pol.','orthogonal pol.'})
ax = axis;
text(ax(1)+0.7*diff(ax(1:2)),ax(3)+0.7*diff(ax(3:4)),{['\tau_1 = ' mnum2str(p(1),1,2) ' ns'], ['\tau_2 = ' mnum2str(p(2),1,2) ' ns']})


