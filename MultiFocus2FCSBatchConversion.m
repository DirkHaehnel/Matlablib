function MultiFocus2FCSBatchConversion(dirname)

fnames = dir([dirname '*.ht3']);

for j=1:length(fnames)
    if ~(exist([dirname fnames(j).name(1:end-3) 'mat'])==2)
        [res, head] = MultiFocus2FCS([dirname fnames(j).name], 4);
        save([dirname fnames(j).name(1:end-4)],'res','head');
    end
end

