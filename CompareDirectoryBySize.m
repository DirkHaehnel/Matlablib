function [r1, r2] = CompareDirectoryBySize(ndir1, ndir2);

dir1 = dir(ndir1);
dir2 = dir(ndir2);

for j=1:length(dir1)
    s1(j) = dir1(j).bytes;
    n1{j} = dir1(j).name;
end
for j=1:length(dir2)
    s2(j) = dir2(j).bytes;
    n2{j} = dir2(j).name;    
end

tmp = setdiff(s1,s2);
r1 = {}; cnt = 1;
for j=1:length(tmp)
    ind = 1:length(s1);
    ind = ind(s1==tmp(j));
    for k=1:length(ind)
        r1{cnt} = n1{ind(k)};
        cnt = cnt+1;
    end
end

tmp = setdiff(s2,s1);
r2 = {}; cnt = 1;
for j=1:length(tmp)
    ind = 1:length(s2);
    ind = ind(s2==tmp(j));
    for k=1:length(ind)
        r2{cnt} = n2{ind(k)};
        cnt = cnt+1;
    end
end

if nargout==0
    for j=1:length(r1)
        disp(r1{j})
    end
end