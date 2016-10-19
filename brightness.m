
names = dir('Point_1.ht3');

ind = ones(numel(names),1);

for n = 1: numel(names)

    FitResults = [];

    load([names(n).name(1:end-4) '.mat'], 'FitResults');

    if ~isempty(FitResults)
        for k = 1:numel(FitResults)
            c = FitResults{1}.c;
            tmp(k,:)  = c(2,:)./sqrt(c(1,:));
        end

        b(n,:) = mean(tmp,1);
 
        l(n) = str2num(names(n).name(7:8));
    else
        ind(n) = 0;
    end

    
end

plot(l(ind==1),b(ind==1,:),'x');
xlabel('Point')
ylabel('brightness / cps)')
