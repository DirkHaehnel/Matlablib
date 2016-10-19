% ScannerCalibration

path = 'D:\Joerg\Doc\Fcs\2Focus\ScannerCalibration\2006-08-17\'

if 0
    names = dir([path 'x*.bmp'])
        
    clear x pos
    for j=1:length(names)
        x(:,:,j) = 256-sum(double(imread([path names(j).name])),3);
        pos(j) = str2num(names(j).name(2:end-4));
    end

    [pos, ind] = sort(pos);
    x = x(:,:,ind);
    x = x(10:end-10,10:end-10,:);

    clear mm
    y = squeeze(sum(x,1));
    for p = 1:10
        tst = [];
        ind = 1:p*10;
        for j=1:size(x,3)-1-p
            for k=1:length(ind)
                tst(k) = sum(y(k+1:end,j).*y(1:end-k,j+p));
            end
            mm(j,p) = mean(ind(tst==max(tst)));
        end
    end
    
    ind = 1:size(y,1);
    for j=1:size(x,3)
        tst = maxl(y(:,j)) & y(:,j)>0.6*mean(y(:,j));
        dd(j) = mean(diff(ind(tst)));
    end
    dd = dd(dd>0.9*mean(dd));    

    for j=1:10
        tmp = mm(1:end-j+1,j);
        xres(j) = mean(tmp(abs(1-tmp/mean(tmp))<0.5))/mean(dd)/j
    end
    
    F = @(p,x,y) sum((y-(p(1)-p(2)*exp(-p(3)*x))).^2);
    px = simplex(F,[0.05 0.05 0.5],[],[],[],[],1:10,xres);
    facx = 200/2^12./px(1);
    
end

if 0
    names = dir([path 'y*.bmp'])

    clear x pos
    for j=1:length(names)
        x(:,:,j) = 256-sum(double(imread([path names(j).name])),3)';
        pos(j) = str2num(names(j).name(2:end-4));
    end

    [pos, ind] = sort(pos);
    x = x(:,:,ind);
    x = x(10:end-10,10:end-10,:);

    clear mm
    y = squeeze(sum(x,1));
    for p = 1:10
        tst = [];
        ind = 1:p*10;
        for j=1:size(x,3)-1-p
            for k=1:length(ind)
                tst(k) = sum(y(k+1:end,j+p).*y(1:end-k,j));
            end
            mm(j,p) = mean(ind(tst==max(tst)));
        end
    end
    
    ind = 1:size(y,1);
    for j=1:size(x,3)
        tst = maxl(y(:,j)) & y(:,j)>0.6*mean(y(:,j));
        dd(j) = mean(diff(ind(tst)));
    end
    dd = dd(dd>0.9*mean(dd));
    
    for j=1:10
        tmp = mm(1:end-j+1,j);
        yres(j) = mean(tmp(abs(1-tmp/mean(tmp))<0.5))/mean(dd)/j
    end

    F = @(p,x,y) sum((y-(p(1)-p(2)*exp(-p(3)*x))).^2);
    py = simplex(F,[0.05 0.05 0.5],[],[],[],[],1:10,yres);
    facy = 200/2^12./py(1);

end

if 1
    path = 'D:\Joerg\Doc\Fcs\2Focus\ScannerCalibration\';
    
    %fname = 'D:\Joerg\Doc\Fcs\2Focus\ScannerCalibration\x_1000x1000_500uspropix.t3r';
    names = dir([path 'x*a.t3r'])
        
    clear p1 p2 zpos
    for j=1:length(names)
        zpos(j) = str2num(names(j).name(3:5));
    end

    [zpos, ind] = sort(zpos);
    names = names(ind);
    
    for j=1:length(names)
       
        if j==1
            [tch, tch, tch, bin] = SCX200Read([path names(j).name]);
            semilogy(bin,tch);
            pos = ginput(2);
        else
            [bin, bin, bin, bin] = SCX200Read([path names(j).name]);
        end
        
        filter = bin>pos(1,1) & bin<pos(2,1);
        [tag, tim1] = SCX200Read([path names(j).name],[filter' filter']);
        tim1 = sum(tim1,3);
        filter = bin<pos(1,1) | bin>pos(2,1);
        [tag, tim2] = SCX200Read([path names(j).name],[filter' filter']);
        tim2 = sum(tim2,3);

        if 0
            close; p1(:,j)=simplex('GridFit',[200 110 20],[],[],[],[],sum(tim1),3)
            close; p2(:,j)=simplex('GridFit',[200 110 20],[],[],[],[],sum(tim2),3)
            for k=1:2
                p1(:,j)=simplex('GridFit',p1(:,j),[],[],[],[],sum(tim1),3)
                p2(:,j)=simplex('GridFit',p2(:,j),[],[],[],[],sum(tim2),3)
            end
        else
            close; p1(:,j)=simplex('GridFit',[200 110 20],[],[],[],[],sum(tim1'),3)
            close; p2(:,j)=simplex('GridFit',[200 110 20],[],[],[],[],sum(tim2'),3)
            for k=1:2
                p1(:,j)=simplex('GridFit',p1(:,j),[],[],[],[],sum(tim1'),3)
                p2(:,j)=simplex('GridFit',p2(:,j),[],[],[],[],sum(tim2'),3)
            end
        end
    end
end

