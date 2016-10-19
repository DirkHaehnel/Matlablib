%[res, head] = SingleFocusRot2FCS('Point_2.ht3');save('Point_2.mat','res','head');

%dirname = 'm:\MTusers\Christoph\Rotation\2011-04-05 Ald Cy5 ulong\'; 
% use only 5 to 15 !!

dirname = 'm:\MTusers\Christoph\Rotation\2011-05-27 Importin Cy5 bis\'; 

if 0
    y = [];
    for i=6:15
        file = ['Point_' num2str(i, '%d')];
        load([dirname file '.mat']);
        [y2, t] = FCSCrossReadRot(res);
        if (isempty(y))
            y = y2;
        else
            y = cat(3, y, y2);
        end
    end
end

nmax = 50;
p2 = zeros(2,nmax); p3 = zeros(3,nmax);
for k=1:nmax
    y = [];
    for kk=1:round(size(y,3)/2)
        file = ['Point_' num2str(round(0.5+4*rand), '%d')];
        %file = ['Point_' num2str(round(4.5+11*rand), '%d')];
        load([dirname file '.mat']);
        [y2, t] = FCSCrossReadRot(res);
        if (isempty(y))
            y = y2;
        else
            y = cat(3, y, y2);
        end
    end
    
    ind = t>1e3;
    for j=1:4
        px = Simplex('ExpFun',1e3,[],[],[],[],t(ind),sum(y(ind,j,:),3),1);
        [err, c] = ExpFun(px, t(ind), sum(y(ind,j,:),3), 1);
        z(:,j) = sum(y(:,j,:),3)./ExpFun(px, c, t)-1;
    end
    
    ind = t<1e3;
    if 1
        p2(:,k) = Simplex('RotoDiffCWFit',[40 500],[],[],[],[], t(ind), z(ind,1:4), [1 1]);
        for kk=1:5
            p2(:,k) = Simplex('RotoDiffCWFit',p2(:,k),[],[],[],[], t(ind), z(ind,1:4), [1 1]);
        end
    end
    if 1
        p3(:,k) = Simplex('RotoDiffCWFit',[40 60 500],[],[],[],[], t(ind), z(ind,1:4), [1 1]);
        for kk=1:5
            p3(:,k) = Simplex('RotoDiffCWFit',p3(:,k),[],[],[],[], t(ind), z(ind,1:4), [1 1]);
        end
    end
    disp(k)
end


% tmp = RotoFit(y(t<2e3,:,:),t(t<2e3),[50 5e2]);
% ExportFigure('FitAll', 'png');


return

    ind = t2>1e3 & t2<5e3;
    for j=1:4
        px = Simplex('ExpFun',1e3,[],[],[],[],t2(ind),sum(y2(ind,j,:),3),1);
        [err, c] = ExpFun(px, t2(ind), sum(y2(ind,j,:),3), 1);
        z(:,j) = sum(y(:,j,:),3)./ExpFun(px, c, t)-1;
    end
    ind = t<1e3;
    if 1
        p2(:,k) = Simplex('RotoDiffCWFit',[40 500],[],[],[],[], t(ind), z(ind,1:4), [1 1]);
        for kk=1:5
            p2(:,k) = Simplex('RotoDiffCWFit',p2(:,k),[],[],[],[], t(ind), z(ind,1:4), [1 1]);
        end
    end
    if 1
        p3(:,k) = Simplex('RotoDiffCWFit',[40 60 500],[],[],[],[], t(ind), z(ind,1:4), [1 1]);
        for kk=1:5
            p3(:,k) = Simplex('RotoDiffCWFit',p3(:,k),[],[],[],[], t(ind), z(ind,1:4), [1 1]);
        end
    end
