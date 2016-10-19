clear all
files = dir('m:\Master-FP\0525\FL*.avi');
mov = aviread(['m:\Master-FP\0525\' files(1).name]);
for j=1:length(mov)
    tst(:,:,j) = mconv2(double(mov(j).cdata),Disk(5));
end
[pthx, pthy, mask] = BrownianMotionTracker(tst);

for j=1:length(pthx{end})-1
    resx{j}=[];
    resy{j}=[];
end
for j=1:length(pthx)
    tmp = pthx{j};
    for k=1:length(tmp)-1
        resx{k} = [resx{k} tmp(k+1:end)-tmp(1:end-k)];
    end
    tmp = pthy{j};
    for k=1:length(tmp)-1
        resy{k} = [resy{k} tmp(k+1:end)-tmp(1:end-k)];
    end
end

for j=1:length(resx)
    sigx(j,1)=mean(resx{j}.^2)-mean(resx{j}).^2;
    sigy(j,1)=mean(resy{j}.^2)-mean(resy{j}).^2;
end

for jf=2:length(files)
    clear tst mov resx resy
    mov = aviread(['m:\Master-FP\0525\' files(jf).name]);
    for j=1:length(mov)
        tst(:,:,j) = mconv2(double(mov(j).cdata),Disk(5));
    end
    [pthx, pthy] = BrownianMotionTracker(tst,[],mask);
    
    for j=1:length(pthx{end})-1
        resx{j}=[];
        resy{j}=[];
    end
    for j=1:length(pthx)
        tmp = pthx{j};
        for k=1:length(tmp)-1
            resx{k} = [resx{k} tmp(k+1:end)-tmp(1:end-k)];
        end
        tmp = pthy{j};
        for k=1:length(tmp)-1
            resy{k} = [resy{k} tmp(k+1:end)-tmp(1:end-k)];
        end
    end
    
    for j=1:length(resx)
        sigx(j,jf)=mean(resx{j}.^2)-mean(resx{j}).^2;
        sigy(j,jf)=mean(resy{j}.^2)-mean(resy{j}).^2;
    end
end