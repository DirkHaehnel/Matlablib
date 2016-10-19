%%
dirname = 'm:\MTusers\Christoph\Dye Photophysics\2011-07-15 Att 647N\';

irfname = 'Reflection long';
fnames = dir([dirname 'Att*.mat']);
int = [10300 1036 107 14270 17270 19800 2013 22640 251 4046 510];

load([dirname irfname]);
irf = t(:,1);
% irf1 = abs(irf1 - mean(irf1(x>x(end)-5)));
% irf2 = abs(irf2 - mean(irf2(x>x(end)-5)));
irf = irf/max(irf);
for j=1:length(fnames)
    load([dirname fnames(j).name]);
    y1(:,j) = t(:,1);
    y2(:,j) = t(:,2);
end
t = x;
[int, srt] = sort(int);
y1 = y1(:,srt);
y2 = y2(:,srt);
%%

%%
mm1 = max(y1);
mm2 = max(y2);
close
%pexc1 = Simplex('PhotophysicsFun',[0 0.1],[0 0],[0 10],[],[],[0 int],[0 mm1],350,640/1.33,150,75e3/60,670/1.33)
pexc1 = Simplex('PhotophysicsSat',0.01,[],[],[],[],[0 int],[0 mm1]);
[~, cexc1] = PhotophysicsSat(pexc1,[0 int],[0 mm1]);
pexc2 = Simplex('PhotophysicsSat',0.01,[],[],[],[],[0 int],[0 mm2]);
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

close; 
photo1 = [0 0 0.001 0.001 0.001 0.001];
irs = irf([3:end-2 1:2]); irs = irs(ind); irs(end-1:end)= irs(end-3);
for j=1:100 
    photo1 = Simplex('PhotophysicsFit',photo1,[0 0 0 0 0 0],[0 0 inf inf inf inf],[],[],int([1 end]),y1(ind,[1 end]),irs/max(irs),pexc1);
end
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

