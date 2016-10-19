function im = DSOMImageRead(name,sze,cnt)

if nargin<3 || isempty(cnt)
    cnt = 2;
end

fnames = dir(name);
tmp = findstr(name,'\');
name = name(1:tmp(end));
for j=1:length(fnames)
    fin = fopen([name fnames(j).name],'r');
    tmp = fread(fin,prod(sze)*cnt,'int16');
    im(:,:,j) = reshape(tmp(2:cnt:end),sze(1),sze(2));
    fclose(fin);
end




