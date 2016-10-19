function [t1, len] = AutodetectTimeGates(tcspcdata, cnum)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tmp = mean(tcspcdata(:,1:min([2 size(tcspcdata,2)])),2);

tmp2 = sort(tmp);
baseline = mean(tmp2(1:ceil(size(tmp2)/5)));
peak = mean(tmp2(ceil(size(tmp2)*0.99):end));

[ind, num] = mCluster(tmp>(baseline+0.15*(peak-baseline)));

tt = 1:size(tcspcdata, 1);
if ((size(tcspcdata, 2) == 4) && (cnum > 2))
    tmp = mean(tcspcdata(:,3:4),2);
    tmp2 = sort(tmp);
    baseline = mean(tmp2(1:ceil(size(tmp2)/5)));
    peak = mean(tmp2(ceil(size(tmp2)*0.99):end));

    [ind2, num2] = mCluster(tmp>(baseline+0.15*(peak-baseline)));
    
    % is ther overlap??
    overlap = zeros(2,1);
    if ((length(num) >= 4) && ((length(num2) >= 2)))
        max_ind = length(num);
        max_ind2 = length(num2);
        for k = 0:3
            overlap(1) = overlap(1) + sum((ind == (max_ind - k)) .* (ind2 == (max_ind2)));
            overlap(2) = overlap(2) + sum((ind == (max_ind - k)) .* (ind2 == (max_ind2 - 1)));
        end
    end

    if ((length(num) < 4) || (overlap(1) == 0) || (overlap(2) == 0))
        for j=1:2
            tst = tt(ind2==length(num2)-j+1);
            ind(tst) = length(num) + 1;
            num = horzcat(num, num2(length(num2)-j+1));
        end
    end
end

for j=1:cnum
    tst = tt(ind==length(num)-j+1); 
    t1(j) = tst(1); 
end
t1 = sort(t1);
t1 = t1-100;
if t1(1)<1 
    t1(1) = 1;
end
t2 = [t1(2:end) max(tt)];
len = min(t2-t1);

end

