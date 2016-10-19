function [xin xhi xlo xof tim] = ReadBenData(name,ind)

fin = fopen(name);
for j=1:ind(1)-1
    if ~feof(fin)
        fgetl(fin);
    end
    j
end

tim = zeros(ind(2),1);
xin = tim;
xhi = tim;
xlo = tim;
xof = tim;

cnt = 1;
tmp = fgetl(fin);
pos = findstr(',',tmp);
s = tmp(pos(2)+1:pos(3)-1);
tim(cnt) = datenum(tmp(pos(1)+1:pos(2)-1),'yyyymmdd') + str2num(s(1:2))/24 + str2num(s(3:4))/1440 + str2num(s(5:6))/86400;
xin(cnt) = str2num(tmp(pos(3)+1:pos(4)-1));
xhi(cnt) = str2num(tmp(pos(4)+1:pos(5)-1));
xlo(cnt) = str2num(tmp(pos(5)+1:pos(6)-1));
xof(cnt) = str2num(tmp(pos(6)+1:end));
for cnt=2:ind(2)
    if ~feof(fin)
        tmp = fgetl(fin);
        s = tmp(pos(2)+1:pos(3)-1);
        tim(cnt) = datenum(tmp(pos(1)+1:pos(2)-1),'yyyymmdd') + str2num(s(1:2))/24 + str2num(s(3:4))/1440 + str2num(s(5:6))/86400;
        xin(cnt) = str2num(tmp(pos(3)+1:pos(4)-1));
        xhi(cnt) = str2num(tmp(pos(4)+1:pos(5)-1));
        xlo(cnt) = str2num(tmp(pos(5)+1:pos(6)-1));
        xof(cnt) = str2num(tmp(pos(6)+1:end));
    else
        break
    end
end

if cnt<ind(2)
    tim(cnt+1:end) = [];
    xin(cnt+1:end) = [];
    xhi(cnt+1:end) = [];
    xlo(cnt+1:end) = [];
    xof(cnt+1:end) = [];
end

fclose(fin);