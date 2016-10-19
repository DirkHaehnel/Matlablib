close all
clear all

mol = [1 2 3 4 6 13 14 15 18 19];
irf = txt2mat('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\SMQauntumYield\Lifetime curves\IRF 10MHz.txt','ReplaceChar',{',.'});

for k=1:length(mol)
    
    dname = ['c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\SMQauntumYield\Lifetime curves\' int2str(mol(k)) '\'];
    fnames = dir([dname '*.txt']);
    
    clear x
    for j=1:length(fnames)
        x(:,:,j)=txt2mat([dname fnames(j).name],'ReplaceChar',{',.'});
        lam(j,k) = str2num(fnames(j).name(1:3));
    end
    
    y = squeeze(x(:,2,:));
    t = squeeze(x(:,1,1));
    
    ind = t>1+t(irf==max(irf)) & t<15+t(irf==max(irf));
    
    for j=1:length(fnames) 
        tau(j,k)=Simplex('ExpFun',1,0,[],[],[],t(ind),y(ind,j),1); 
        for jj=1:3
            tau(j,k)=Simplex('ExpFun',tau(j,k),0,[],[],[],t(ind),y(ind,j),1); 
        end
    end
    
    plot(lam.^2./lam,tau)
    
end