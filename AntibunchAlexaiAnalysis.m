flag = 1;

if 1
   tmp = conv(res.rate,ones(1,3))/3;
   ind = abs((tmp(2:end-1)-res.rate)./tmp(2:end-1))<=0.05;
   res.auto = res.auto(:,:,ind);
   res.rate = res.rate(ind);
   res.time = res.time(ind);
   anti.auto = anti.auto(:,:,ind);
end

mm = 10;
ng = ceil(size(res.auto,3)/mm);

close all

% ind = res.autotime>0; 
% ind(1) = false;
% clear c
% p0 = simplex('rigler',[3e-2 5e-4 1e5 1e6],[],[],[],[],res.autotime(ind),mean(sum(res.auto(ind,:,:),3),2));
% for k=1:2
%     p0 = simplex('rigler',p0,[],[],[],[],res.autotime(ind),mean(sum(res.auto(ind,:,:),3),2));
% end
% for j=1:mm
%     p = simplex('rigler',p0,[],[],[],[],res.autotime(ind),mean(sum(res.auto(ind,:,(j-1)*ng+1:min([j*ng,end])),3),2));
%     for k=1:2
%         p = simplex('rigler',p,[],[],[],[],res.autotime(ind),mean(sum(res.auto(ind,:,(j-1)*ng+1:min([j*ng,end])),3),2));
%     end
%     [err,c(:,j)] = rigler(p,res.autotime(ind),mean(sum(res.auto(ind,:,(j-1)*ng+1:min([j*ng,end])),3),2));
% end
% lam = c(2,:)./sum(c(2:end,:));

ind = res.autotime<=1e-5; 
ind(1) = false;
clear c off
tst0 = [0.03 0.003 3.5e3 4e4 2e5];
for j=1:mm
    tst = tst0
    for k=1:3
        tst = simplex('Rigler',tst,0*tst,[],[],[],res.autotime(2:end),sum(sum(res.auto(2:end,:,(j-1)*ng+1:min([j*ng,end])),3),2),[],1)
    end
    [err,c(:,j)]=Rigler(tst,res.autotime(2:end),sum(sum(res.auto(2:end,:,(j-1)*ng+1:min([j*ng,end])),3),2),[],1);
    off(j) = mean(sum(sum(res.auto(end-10:end,:,(j-1)*ng+1:min([j*ng,end])),3),2));
end

ind = anti.autotime<=37.5;%125; 
clear cc c0
pp0 = simplex('AntiBunch',[0 1],[],[],[],[],anti.autotime(ind),sum(anti.auto(ind,:,:),3),25,[],1);
for k=1:2
    pp0 = simplex('AntiBunch',pp0,[],[],[],[],anti.autotime(ind),sum(anti.auto(ind,:,:),3),25,[],1);
end
for j=1:mm
    pp = simplex('AntiBunch',pp0,[],[],[],[],anti.autotime(ind),sum(anti.auto(ind,:,(j-1)*ng+1:min([j*ng,end])),3),25);
    for k=1:2
        pp = simplex('AntiBunch',pp,[],[],[],[],anti.autotime(ind),sum(anti.auto(ind,:,(j-1)*ng+1:min([j*ng,end])),3),25);
    end
    [err,cc(:,j),c0(:,j),z,z,z] = AntiBunch(pp,anti.autotime(ind),sum(anti.auto(ind,:,(j-1)*ng+1:min([j*ng,end])),3),25);
end

%lam = c(1,:)./sum(c);
% h = (c(1,:)-off)./sum(c).*sum(cc(1:2,:));
% delta = (lam-1).*sum(cc(1:2,:)) + cc(3,:);
% h = (c(1,:)-off)./sum(c).*cc(2,:);
% lam = (c(1,:)-off)./(sum(c)-off);
% h = lam.*(cc(2,:)-c0(2,:)*60/57);
% delta = (c(1,:)./sum(c)-1).*cc(2,:) + cc(3,:);
% nemit = h./delta;
% npart = c0(2,:)./h;

lam = c(2,:)./sum(c(2:end,:))
h = lam.*(cc(2,:)-c0(2,:)*60/57)
delta = (sum(c(1:2,:))./sum(c)-1).*cc(2,:) + cc(3,:)
nemit = h./delta


return

for cnt=1:41 eval(['load Data' mint2str(cnt,2)]); AntibunchMultiAnalysis; final(cnt).nemit = nemit; final(cnt).h = h; final(cnt).delta = delta; final(cnt).c = c; final(cnt).cc = cc; final(cnt).off = off; final(cnt).lam = lam; save AntibunchResults names final; end

