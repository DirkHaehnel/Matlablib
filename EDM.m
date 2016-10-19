function tt = EDM(t);

% tt = EDM(t) calculates the Euclidian Distance Map of array t
% (C) Joerg Enderlein, 2001

[m,n] = size(t);
tt = (t>0);
tmpold = 0;
tmp = sum(sum(tt));
ind = tt>0;
while ~(tmpold==tmp)
   tmpold = tmp;
   tmp = min(min(min(tt(:,[1 1:end-1]),tt([1 1:end-1],:)),tt(:,[2:end end])),tt([2:end end],:))+2;
   tmpd = min(min(min(tt([1 1:end-1],[1 1:end-1]),tt([1 1:end-1],[2:end end])),...
      tt([2:end end],[2:end end])),tt([2:end end],[1 1:end-1]))+3;
   tmp = min(tmp,tmpd);
   tt = tmp.*ind;
   tmp = sum(sum(tt));
end

