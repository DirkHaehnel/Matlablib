%%
dirname = 'm:\MTusers\Christoph\Dye Photophysics\2011-06-14 Alx 647\';

irfname = 'Milk 8500';
fnames = dir([dirname 'Alx*.mat']);
int = [10600 1100 13500 2850 310 52 6350];

load([dirname irfname]);
irf1 = t(:,1);
irf2 = t(:,2);
% irf1 = abs(irf1 - mean(irf1(x>x(end)-5)));
% irf2 = abs(irf2 - mean(irf2(x>x(end)-5)));
irf1 = irf1/max(irf1);
irf2 = irf2/max(irf2);
for j=1:length(fnames)
    load([dirname fnames(j).name]);
    y1(:,j) = t(:,1);
    y2(:,j) = t(:,2);
end
[int, srt] = sort(int);
y1 = y1(:,srt);
y2 = y2(:,srt);
%%

%%
mm1 = max(y1);
mm2 = max(y2);
close
%pexc1 = Simplex('PhotophysicsFun',[0 0.1],[0 0],[0 10],[],[],[0 int],[0 mm1],350,640/1.33,150,75e3/60,670/1.33)
pexc1 = Simplex('PhotophysicsSat',0.1,[],[],[],[],[0 int],[0 mm1]);
[~, cexc1] = PhotophysicsSat(pexc1,[0 int],[0 mm1]);
pexc2 = Simplex('PhotophysicsSat',0.1,[],[],[],[],[0 int],[0 mm2]);
[~, cexc2] = PhotophysicsSat(pexc2,[0 int],[0 mm2]);

t = (1:size(y1,1))'*0.0125;
ind = sum(y1,2)>0.5*mean(sum(y1,2));
t = t - min(t(ind));
%ind = ind & t>5;
clear lam lam1 lam2
for j=1:size(y1,2)
    if 0
        lam(j) = Simplex('ExpFun',[5],[0],[inf],[],[],t(ind),y1(ind,j)+y2(ind,j),1);
        lam(j) = Simplex('ExpFun',lam(j),[0],[inf],[],[],t(ind),y1(ind,j)+y2(ind,j),1);
        lam1(j) = Simplex('ExpFun',[5],[0],[inf],[],[],t(ind),y1(ind,j),1);
        lam1(j) = Simplex('ExpFun',lam1(j),[0],[inf],[],[],t(ind),y1(ind,j),1);
        lam2(j) = Simplex('ExpFun',[5],[0],[inf],[],[],t(ind),y2(ind,j),1);
        lam2(j) = Simplex('ExpFun',lam2(j),[0],[inf],[],[],t(ind),y2(ind,j),1);
    else
        lam(:,j) = Simplex('ExpFun',[1 5],[0 0],[inf inf],[],[],t(ind),y1(ind,j)+y2(ind,j),1);
        lam(:,j) = Simplex('ExpFun',lam(:,j),[0 0],[inf inf],[],[],t(ind),y1(ind,j)+y2(ind,j),1);
        lam1(:,j) = Simplex('ExpFun',[1 5],[0 0],[inf inf],[],[],t(ind),y1(ind,j),1);
        lam1(:,j) = Simplex('ExpFun',lam1(:,j),[0 0],[inf inf],[],[],t(ind),y1(ind,j),1);
        lam2(:,j) = Simplex('ExpFun',[1 5],[0 0],[inf inf],[],[],t(ind),y2(ind,j),1);
        lam2(:,j) = Simplex('ExpFun',lam2(:,j),[0 0],[inf inf],[],[],t(ind),y2(ind,j),1);
    end
end

photo1 = Simplex('PhotophysicsFit',[0.03 0.003],[0 0],[],[],[],int,y1(ind,:),irf1(ind),pexc1);
photo1 = Simplex('PhotophysicsFit',photo1,[0 0],[],[],[],int,y1(ind,:),irf1(ind),pexc1);
p1 = Simplex('ExpFun',[500 2e3],[0 0],[],[],[],1:4.5e3,tcspc(1.5e3+1:6e3,1),1);
[err, c1] = ExpFun(p1,1:4.5e3,tcspc(1.5e3+1:6e3,1));
tauf1 = sum(c1(2:3).*p1)/sum(c1(2:3))*2e-3;

clf
photo2 = Simplex('PhotophysicsFit',photo1,[0 0],[],[],[],int,y2(ind,:),irf2(ind),pexc2);
photo2 = Simplex('PhotophysicsFit',photo2,[0 0],[],[],[],int,y2(ind,:),irf2(ind),pexc2);
p2 = Simplex('ExpFun',[500 2e3],[0 0],[],[],[],1:4.5e3,tcspc(1.5e3+1:6e3,2),1);
[err, c2] = ExpFun(p2,1:4.5e3,tcspc(1.5e3+1:6e3,2));
tauf2 = sum(c2(2:3).*p2)/sum(c2(2:3))*2e-3;

[[1/(photo1(1)/12.5e-9) 1/(tauf1*photo1(2)/12.5/12.5e-9)]*1e6;
[1/(photo2(1)/12.5e-9) 1/(tauf2*photo2(2)/12.5/12.5e-9)]*1e6]
%%

