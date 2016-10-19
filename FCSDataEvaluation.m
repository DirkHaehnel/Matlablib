p1 = 0.02;
w1 = 350; 
a1 = 150;
clear dc v conc w0 a0 triplet c velo err
global pd
for j=2:length(fnames)
    load([dirname fnames(j).name])
    [data.y,data.t] = FCSCrossRead(res,[0 1]); % converting raw correlation functions in what we need
    
    if 0
        tmp = mean(squeeze(data.y(:,1,:)),2);
        covco = zeros(1,size(data.y,3));
        for k=1:size(data.y,3)
            coef = tmp\squeeze(data.y(:,1,k));
            covco(k) = sum(squeeze(data.y(:,1,k)).*(coef*tmp))/sqrt(sum(squeeze(data.y(:,1,k)).^2).*sum((coef*tmp).^2));
        end
        data.y(:,:,covco<0.999)=[];
    end
    
    pd = [p1 1]; % setting initial diffusion and triplet
    [dc{j} v{j} conc{j} w0{j} a0{j} triplet{j} c{j} velo{j} err{j}] = FCSFit(data,[w1 a1],1,1,[75e3/60 [470 520]/1.33 384],[[0 0];[inf 50]]);
    eval(['print -dpng -r300 ' dirname fnames(j).name(1:end-4)])
    eval(['save ' dirname 'Results.mat dc v conc w0 a0 triplet c velo err'])
end
