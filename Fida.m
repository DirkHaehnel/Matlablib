if 1
    cnt = 1;
    bin = 0:100;
    [y, tcspc, chan, markers, num, overcount, head] = pt3v2read('D:\5nM_200mMNaCl.pt3', [1 1e6]);
    cnt = cnt + num;
    
    [auto, autotime] = tttr2xfcs(y,ones(size(y)),20,10);
    while num>0
        [y, tcspc, chan, markers, num, overcount, head] = pt3v2read('D:\5nM_200mMNaCl.pt3', [cnt 1e6]);
        cnt = cnt + num;
        auto = auto + tttr2xfcs(y,ones(size(y)),20,10);
        semilogx(autotime, auto); drawnow
    end
end

if 0
    names = dir('d:\*.pt3');

    for j=1:length(names)
        fida = zeros(101,1);
        num = 1;
        cnt = 1;
        bin = 0:100;
        while num>0
            [y, tcspc, chan, markers, num, overcount, head] = pt3v2read(['d:\' names(j).name], [cnt 1e6]);
            cnt = cnt + num;
            fida = fida + mHist(tttr2bin(y,5e3),bin);
            semilogy(bin, fida); drawnow
        end
        save(names(j).name(1:end-4),'fida','bin');
    end
end