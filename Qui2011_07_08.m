dirname = 'm:\MTusers\Qui\110708_PET-HTH\';
filenames = dir([dirname '*.ht3']);

if 0 % compute correlations
    for j=1:length(filenames)
        [res, head] = TwoFocus2FCS([dirname filenames(j).name], [1 1e2]);
        eval(['save ' dirname '' filenames(j).name(1:end-4) ''' res head']);
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
    
