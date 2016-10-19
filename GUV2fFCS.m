% GUV-2fFcs

dirs ={'D:\Joerg\Doc\Fcs\2Focus\GUV\Glucose\2006-05-08\',...
    'D:\Joerg\Doc\Fcs\2Focus\GUV\Succrose\2006-04-26\',...
    'D:\Joerg\Doc\Fcs\2Focus\GUV\Succrose+Salt\2006-05-02\',...
    'D:\Joerg\Doc\Fcs\2Focus\GUV\Water\2006-04-25\'};

for j=1:length(dirs)
    names = dir([dirs{j} '*parts.mat']);
    for k=1:length(names)    
        [j k]
        eval(['load ' dirs{j} names(k).name ' parameters']);
        tt{j}(k) = str2num(parameters.Temperature);
    end
end



return
global pd
for j=1:length(dirs)
    names = dir([dirs{j} '*parts.mat']);
    clear dc v conc w0 a0
    for k=1:length(names)
        pd = 0.1;
        %[dc(:,k) v(:,k) conc(:,k) w0(:,k)] = FCSFit([dirs{j} names(k).name],[400 0],[],1,[[0 0];[inf 0]]);
        [dc(:,k) v(:,k) conc(:,k) w0(:,k) a0(:,k)] = FCSFit([dirs{j} names(k).name],[400 10],[],1);
    end
    eval(['save ' dirs{j} 'fitres1 dc v conc w0 a0 names'])
end

