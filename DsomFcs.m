% DsomFcs analysis

clear 

% names = {'c:\Users\Public\Transfer\DSOM_FCS\20to200.pt3',...
%     'c:\Users\Public\Transfer\DSOM_FCS\20to200.pt3',...
%     'c:\Users\Public\Transfer\DSOM_FCS\20to200.pt3'};
% 
% if 1
%     for j=1:length(names)
%         [res{j} head{j}] = TwoFocus2FCS(names{j});
%         save 'c:\Users\Public\Transfer\DSOM_FCS\results' res head;
%     end
% end

load c:\Joerg\Doc\Fcs\DSOMFCS\081023_DSOM_FCS\Data.mat

sv = 0:0.01:2;
for j=1:length(fnames)

    [y,t] = FCSCrossRead(res(j));
    y(:,:,end)=[];
    rate = res(j).rate(1:end-1,:)';
    if mean(y(end,1,:),3)>mean(y(end,2,:),3)
        y = y(:,[2 1 3 4],:);
        rate = rate([2 1 4 3],:);
    end
    y(:,1,:) = squeeze(y(:,1,:))./(ones(length(t),1)*(rate(1,:).*rate(3,:)));
    y(:,2,:) = squeeze(y(:,2,:))./(ones(length(t),1)*(rate(2,:).*rate(4,:)));
    y(:,3,:) = squeeze(y(:,3,:))./(ones(length(t),1)*(rate(1,:).*rate(4,:)));
    y(:,4,:) = squeeze(y(:,4,:))./(ones(length(t),1)*(rate(2,:).*rate(3,:)));
    y = mean(y,3);

    y(:,3) = sum(y(:,3:4),2);
    y(:,4) = [];

    p1(:,j) = [1e-3 1e-2]';
    close; for k=1:3 p1(:,j) = Simplex('Rigler',p1(:,j),[0 0],[],[],[],t,y(:,1),[],1); end
    [err,c,z1(:,j)] = Rigler(p1(:,j),t,sum(y(:,1),2));
    p2(:,j) = p1(:,j);
    close; for k=1:3 p2(:,j) = Simplex('Rigler',p2(:,j),[0 0],[],[],[],t,y(:,2),[],1); end
    [err,c,z2(:,j)] = Rigler(p2(:,j),t,sum(y(:,2),2));
    p3(:,j) = p1(:,j);
    close; for k=1:3 p3(:,j) = Simplex('Rigler',p3(:,j),[0 0],[],[],[],t,y(:,3),[],1); end
    [err,c,z3(:,j)] = Rigler(p3(:,j),t,sum(y(:,3),2));

    %z1(:,j)=y(:,1);z2(:,j)=y(:,2);z3(:,j)=y(:,3);
    for k=1:length(sv)
        s=sv(k);
        zz=z1(:,j)+s^2*z2(:,j)-s*z3(:,j);
        semilogx(t,zz); drawnow
        tau(k,j) = interp1(zz,t,(1-s)^2+0.5*(max(zz)-(1-s)^2),'cubic');
        pn(k,j) = zz(end)/(zz(1)-zz(end));
        int(k,j) = zz(end);
    end
    
%    save 'c:\Joerg\Doc\Fcs\DSOMFCS\DSOMFcsData' sv z1 z2 z3 sv tau pn int
end

return 

plot(sv,1e6*tau./(ones(size(sv'))*tau(1,:))*tau(1,1))
legend({'20 vs 200 \muW','20 vs 300 \muW','20 vs 400 \muW'},2)
xlabel('subtracting factor'); ylabel('diffusion time (\mus)')

plot(sv,tau./(ones(size(sv'))*tau(1,:)))
legend({'20 vs 200 \muW','20 vs 300 \muW','20 vs 400 \muW'},2)
xlabel('subtracting factor'); ylabel('rel. diffusion time')

[ind,ind] = min(tau);
tmp = pn./(ones(size(sv'))*pn(1,:));
plot(sv,tmp,sv(ind(1)),tmp(ind(1),1),'or',sv(ind(2)),tmp(ind(2),2),'ob',sv(ind(3)),tmp(ind(3),3),'og')
axis([0 2 0 5])
legend({'20 vs 200 \muW','20 vs 300 \muW','20 vs 400 \muW'},2)
xlabel('subtracting factor'); ylabel('rel. particle number')

