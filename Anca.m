% Program for FLCS analysis of Anca's data

fluid = {'07dopc','06dopc'};
solid = {'06lo','07lo','08lo'};
mix = {'14LO-LD','18LO-LD'};

if 1 % data conversion
    names = solid;
    for k = 2:2%length(names)
        tmp = dir(['D:\Daten\AncaMargineanu\' names{k} '*.spc']);
        x = []; y = [];
        for j=1:length(tmp)
            q = spcreader(['D:\Daten\AncaMargineanu\' tmp(j).name]);
            y = [y; q(:,2)];
            x = [x; q(:,2)];
            eval(['save Anca' names{k} ' x y'])
        end
    end
end

if 1 % lifetime analysis
    cnt = 1;
    names = fluid;
    bin = 0:255;
    for j=1:length(names)
        eval(['load Anca' names{j} ' x']);
        h(:,cnt) = mhist(x,bin);
        cnt = cnt+1;
    end
    names = solid;
    for j=1:length(names)
        eval(['load Anca' names{j} ' x']);
        h(:,cnt) = mhist(x,bin);
        cnt = cnt+1;
    end
    names = mix;
    for j=1:length(names)
        eval(['load Anca' names{j} ' x']);
        h(:,cnt) = mhist(x,bin);
        cnt = cnt+1;
    end
    h = h(end:-1:1,:);
    t = repmat(bin',[1 size(h,2)]);
    t(:,3) = 0.14108828*t(:,3); 
    t(:,[1:2 4:end]) = 0.13020833*t(:,[1:2 4:end]);
    save AncaLifetimes t h
    tmin = 50; tmax = 240;
    ind = bin>=tmin & bin<=tmax;
    for j=1:size(h,2)
        p(:,j) = simplex('ExpFun',[1 2],[0 0],[],[],[],t(ind,j),h(ind,j),1); 
    end
    for j=1:size(h,2)
        p1(:,j) = simplex('ExpFun',2,0,[],[],[],t(ind,j),h(ind,j),1); 
    end
    save AncaLifetimes p -append
    
    semilogy(t,(h-ones(size(h,1),1)*mean(h(1:10,:)))./(ones(size(h,1),1)*(h(50,:)-mean(h(1:10,:)))))
    legend({'07dopc','06dopc','06lo','07lo','08lo','14LO-LD','18LO-LD'},4)
    axis([0 30 1e-5 10])
    xlabel('time [ns]'); ylabel('rel. intensity')
end

