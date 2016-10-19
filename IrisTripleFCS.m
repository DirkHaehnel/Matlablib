av = 100e3/60;
lam = [640 670]/1.33;
dist = 403;

load c:\Joerg\Doc\Fcs\2Focus\Recoverin\wtmRecX647_10mgmlDOPC_250C_P18.mat

data.y(:,1,:) = sqrt(res.auto(:,1,3,:).*res.auto(:,3,1,:));
data.y(:,2,:) = sqrt(res.auto(:,2,4,:).*res.auto(:,4,2,:));
data.y(:,3,:) = sqrt(res.auto(:,1,4,:).*res.auto(:,3,2,:));
data.y(:,4,:) = sqrt(res.auto(:,2,3,:).*res.auto(:,4,1,:));
data.t = res.autotime;

global pd
pd = [0.1 0.001 10 20];
[dc v conc w0 a0 triplet c1] = FCSFit(data,[350 150],length(pd)-2,[],[av lam dist]);
[dc v conc w0 a0 triplet c1] = FCSFit(data,[w0 a0],length(pd)-2,[],[av lam dist]);

[zd zd] = TripleFCS(res); zd = sum(sum(zd,3),2);

bounds = [pd(1:2) 0 0; pd(1.2) inf inf];
para = Simplex('Gauss3Fcs',[w0 a0],[],[],[],[],av,lam,res.autotime*1e6,zd,2,bounds);

zm = Gauss3FCS([w0 a0],av,lam,res.autotime*1e6/pd(1));
zm = [zm Gauss3FCS([w0 a0],av,lam,res.autotime*1e6/pd(2))];
zm = [exp(-res.autotime*1e6*(1./triplet')) zm];
c2 = lsqnonneg(zm,zd);

semilogx(res.autotime,zd,'o',res.autotime,zm*c2)
semilogx(res.autotime,zd,'o',res.autotime,cumsum(zm(:,[3 4 1 2]).*(ones(length(res.autotime),1)*c2([3 4 1 2])'),2))

cc = mean(c1(end-1:end,:),2).^3./c2(end-1:end).^2;