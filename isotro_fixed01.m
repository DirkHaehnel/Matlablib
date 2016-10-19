close all
clear
timeunit=0.1

[filename, filepath] = uigetfile('*.mat');
tit=findstr(filename, '.');tit=filename(1:tit(end)-1);
filename = [filepath filename];
outfile = findstr(filename,'.');outfile = filename(1:outfile(end)-1);
eval(['load ' filename ' res'])


bin = 0:0.5:180;
for mol=1:1%54
    
    kv=1:60; 
    %eval(['load ' outfile 'M' mint2str(mol,3) ' res'])
   
    clear hh pp
    bin = 0:0.5:50;
    theta = res(:,4)/180*pi;
    phi = unwrap(res(:,5)/180*2*pi)/2;
    theta = mod(round(res(:,5)/180-phi/pi),2)+(-1).^round(res(:,5)/180-phi/pi).*theta/pi;
    t = 1:size(res,1);
    t = t(~(res([1:end-1],1)+1==res([2:end],1)));
    t = [0 t size(res,1)];
    for j=1:length(t)-1
        pp(j).theta = theta(t(j)+1:t(j+1));
        pp(j).phi = phi(t(j)+1:t(j+1));
    end
    
    j = 1;
    tmp = [sin(pp(j).theta).*cos(pp(j).phi), sin(pp(j).theta).*sin(pp(j).phi), cos(pp(j).theta)];
    for k=1:max(kv)
        if k<size(tmp,1)
            hh(:,k) = mhist(acos(sum(tmp(1:end-k,:).*tmp(k+1:end,:),2))/pi*180,bin);
            mm(k,mol) = sum(sum(tmp(1:end-k,:).*tmp(k+1:end,:),2));
            cnt(k) = size(tmp,1)-k;
        else
            hh(:,k) = zeros(length(bin),1);
            mm(k,mol) = 0;
            cnt(k) = 0;
        end            
    end
    for j=2:length(pp)
        tmp = [sin(pp(j).theta).*cos(pp(j).phi), sin(pp(j).theta).*sin(pp(j).phi), cos(pp(j).theta)];
        for k=1:max(kv)
            if k<size(tmp,1)
                hh(:,k) = hh(:,k) + mhist(acos(sum(tmp(1:end-k,:).*tmp(k+1:end,:),2))/pi*180,bin);
                mm(k,mol) = mm(k,mol) + sum(sum(tmp(1:end-k,:).*tmp(k+1:end,:),2));
                cnt(k) = cnt(k) + size(tmp,1)-k;
            end
        end
    end
    mm(:,mol) = mm(:,mol)./cnt';
    for ii=1:size(hh,2)
        hh(:,ii)=hh(:,ii)./cnt(ii)';
    end
    z(:,mol) = (cos(bin/180*pi)*hh./(sum(hh)+(sum(hh)==0)))';
    z2(:,mol) = ((cos(bin/180*pi).^2)*hh./(sum(hh)+(sum(hh)==0)))';
   
    tph=['M' mint2str(mol,3)];
    
    %<cos>
    figure(1);plot(z(:,mol), 'o'); title(tph, 'fontsize',15);
    coef(:,mol) = simplex('ExpPur',[1e0 1e-2],[0 0],[],[],[],[0 kv(1:30)],[1; z(z(1:30,mol)>0,mol)]);
    for s=1:3
        coef(:,mol) = simplex('ExpPur',coef(:,mol),[0 0],[],[],[],[0 kv(1:30)],[1; z(z(1:30,mol)>0,mol)]);
    end
        [err(mol), weight(:,mol), tmp, zz(:,mol)] = ExpPur(coef(:,mol),[0 kv(1:30)],[1; z(z(1:30,mol)>0,mol)]);
    end
     
    Drot(:,mol)=coef(:,mol)/2/timeunit
    WT1(:,mol)=weight(:,mol)
    
    
    %<cos>^2
     figure(2);plot(z2(:,mol), 'o'); title(tph, 'fontsize',15); 
    coef2(:,mol) = simplex('ExpPur',[1e0 1e-2 1e-4],[0 0 0],[],[],[],[0 kv],[1; z2(z2(:,mol)>0,mol)]);
    for s=1:3
        coef2(:,mol) = simplex('ExpPur',coef2(:,mol),[0 0 0],[],[],[],[0 kv],[1; z2(z2(:,mol)>0,mol)]);
    end
        [err2(mol), weight2(:,mol), tmp2, zz2(:,mol)] = ExpPur(coef2(:,mol),[0 kv],[1; z2(z2(:,mol)>0,mol)]);
    end
   Drot2(:,mol)=coef2(:,mol)/2/timeunit
    WT2(:,mol)=weight2(:,mol)

