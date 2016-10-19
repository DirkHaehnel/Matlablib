% FrankFitting

close all
clear all

files = dir('c:\Joerg\Doc\Meixner\Frank\Härtekinetik\*.txt');

for j=1:length(files)
    x = load(['c:\Joerg\Doc\Meixner\Frank\Härtekinetik\' files(j).name]);
    ind = 1:size(x,1);
    clear y
    for k=1:size(x,2)-1
        pos = round(mean(ind(x(:,k+1)==max(x(:,k+1)))));
        y(k) = x(pos,1);
    end
    len(j) = y(end);
    t = (1:size(x,2)-1)*0.5;
    tst=conv(y,ones(1,5))/5;
    tst = abs(y-tst(3:end-2))./y>0.02;
    tst(1) = false; tst(end) = false;
    p(j) = Simplex('ExpFun',30,0,[],[],[],t(~tst),y(~tst),1); 
    err(j) = ExpFun(p(j),t(~tst),y(~tst));
end
    
[len, ind] = sort(len);
p = p(ind);

