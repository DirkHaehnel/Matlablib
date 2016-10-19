clear xc yc
data(:,1)=round(data(:,1)./10)+5;
data(:,2)=round(data(:,2)./10)+5;
for i=1:500;
    ind=data(:,3)==(i);
    xc{i}=data(ind,1)';
    yc{i}=data(ind,2)';
end

    