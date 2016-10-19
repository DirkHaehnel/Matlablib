dirname = 'm:\MTusers\Qui\110719_Atto in GdnHCl\';
filenames = dir([dirname '*.ht3']);

if 1 % compute correlations
    for j=2:length(filenames)
        [res, head] = SingleFocus2FCS([dirname filenames(j).name], [1 1e2]);
        eval(['save ''' dirname '' filenames(j).name(1:end-4) ''' res head']);
    end
end

dirname = 'm:\MTusers\Qui\110721_PET_HTH_Fxa\';
filenames = dir([dirname '*.ht3']);

if 1 % compute correlations
    for j=1:length(filenames)
        [res, head] = SingleFocus2FCS([dirname filenames(j).name], [1 1e2]);
        eval(['save ''' dirname '' filenames(j).name(1:end-4) ''' res head']);
    end
end

dirname = 'm:\MTusers\Qui\110728_PET_Test\';
filenames = dir([dirname '*.ht3']);

if 1 % compute correlations
    for j=1:length(filenames)
        [res, head] = SingleFocus2FCS([dirname filenames(j).name], [1 1e2]);
        eval(['save ''' dirname '' filenames(j).name(1:end-4) ''' res head']);
    end
end

if 0
    for j=1:length(filenames)
        eval([dirname filenames(j).name(1:end-4)]);
        [tmp t] = FCSCrossRead(res);
        z(:,j) = sum(sum(tmp,3),2);
    end
end

close; 
for j=1:length(filenames)
    p = Simplex('ExpFun',1e-6,0,[],[],[],t(ind),y(ind,j)/mean(y(end-10:end,j))-1,1,1);
    subplot(4,1,1:3)
    title(['\tau = ' mnum2str(1e6*p,1,2) ' \mus']);
    eval(['print -dpng -r300 ' dirname filenames(j).name(1:end-4)]);
end
    
