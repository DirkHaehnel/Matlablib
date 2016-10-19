function z = FilterFCS(t,y,nn,th)

if nargin<3 || isempty(nn)
    nn = 30;
end
if nargin<4 || isempty(th)
    th = 2;
end

z = zeros(size(y,1),size(y,2),size(y,3));
for j=1:size(y,2)
    tmp = squeeze(y(:,j,:));
    mm = mean(tmp);
    ind = mm<th*median(mm);
    mm = (ones(size(tmp,1),1)*mm);
    tmp = tmp - mm;
    sc = sqrt(sum(tmp.^2));
    ind = ind & sc<th*median(sc);
    sc = (ones(size(tmp,1),1)*sc);
    tmp = tmp./sc;
    for jj=1:size(tmp,2)
        c(:,jj) = polyfit(t(end-nn:end),tmp(end-nn:end,jj),1);
    end
    ind = ind & abs(c(1,:))<th*mean(abs(c(1,:)));
    tmp = tmp.*sc + mm;
    z(:,j,ind) = tmp(:,ind);
end

