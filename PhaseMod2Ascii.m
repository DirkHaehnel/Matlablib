function head = PhaseMod2Ascii(fname)

photons = 5e5; % number of photons processed one at a time

fout = fopen([fname(1:end-3) 'dat'],'w');

cnt = 0;
tst = 1;
[head, head, head, head] = Spcread_4096(fname,[0 0]);
Timeunit = head.GlobClock*1e-9;
while tst
    [y, flag, num, num, num] = Spcread_4096(fname, [cnt+1 photons]);
    fprintf(fout, '%u %u\n', [y flag]')

    cnt = cnt + num;
    if num<photons
        tst = 0;
    end
end

fclose(fout);

