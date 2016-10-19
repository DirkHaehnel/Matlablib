fname ='c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\decays Alexey\IRF.txt';
irf = load(fname);

dname = 'c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\decays Alexey\R6G in water\';
fnames = dir([dname '*.txt']);
for j=1:length(fnames)
    x1(:,:,j) = load([dname fnames(j).name]);
end
for j=1:size(x1,3)-2 
    lam1(j)=str2num(fnames(j).name(1:3)); 
end

if 1
    ind = x1(:,1,1)>1.8 & x1(:,1,1)<12;
    close; 
    for j=1:size(x1,3) 
        p1(:,j)=Simplex('ExpFun',1,0,[],[],[],x1(ind,1,j),x1(ind,2,j),2); 
    end
    plot(lam1,p1(1:end-2),lam1,p1(end-1)*ones(size(lam1)),':');
end
    
dname = 'c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\decays Alexey\R6G in glycerol\';
fnames = dir([dname '*.txt']);
for j=1:length(fnames)
    x2(:,:,j) = load([dname fnames(j).name]);
end
for j=1:size(x2,3)-1
    lam2(j)=str2num(fnames(j).name(1:3)); 
end
    
if 1
    ind = x2(:,1,1)>1.8 & x2(:,1,1)<12;
    close; 
    for j=1:size(x2,3) 
        p2(:,j)=Simplex('ExpFun',1,0,[],[],[],x2(ind,1,j),x2(ind,2,j),2); 
    end
    plot(lam2,p2(1:end-1),lam2,p2(end)*ones(size(lam2)),':');
end

plot(lam1,p1(1:end-2),lam1,p1(end-1)*ones(size(lam1)),':',lam2,p2(1:end-1),lam2,p2(end)*ones(size(lam2)),':');
xlabel('transmission max. (nm)')
ylabel('lifetime (ns)')
legend({'water','free water','glycerol','free glycerol'},2)
