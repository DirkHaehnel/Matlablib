function tim = FLIM2Mean(flim,tcspc)

if nargin>1 || ~isempty(tcspc)
    [pos,pos] = max(tcspc);
    %bck = mean(tcspc(2:max([2 pos-10])));
else
    pos = 1; bck = 0;
end
for j=1:size(flim,1)
    for k=1:size(flim,2)
        len = sum(flim{j,k}>=pos-1);
        if len>0
            tim(j,k) = mean(flim{j,k}(flim{j,k}>=pos-1)-pos+1);
        else
            tim(j,k) = 0;
        end
    end
end
