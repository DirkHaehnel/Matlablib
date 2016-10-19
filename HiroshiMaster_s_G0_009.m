% HiroshiMaster

clear all

fname1 = {'D:\Daten\Johan\s_G0_009\s_g0_009_1-500.TIF', 'D:\Daten\Johan\s_G0_009\s_g0_009_501-1000.TIF', ...
    'D:\Daten\Johan\s_G0_009\s_g0_009_1001-1500.TIF', 'D:\Daten\Johan\s_G0_009\s_g0_009_1501-1700.TIF'};
outfile = 'D:\Daten\Johan\s_G0_009\s_g0_009';

if 0
    filename = fname1;
    imnum = length(imfinfo(filename));
    % localize regions of interest
    im = double(imread(filename,1));
    for j=2:imnum
        im = im + double(imread(filename,j));
    end
    mim(im);
    pos = ginput;
    pos = round(pos);
end

pos = [
   228    26
   168    39
   132    61
   137    91
   342    69
   469    95
   257    86
   290    93
   296   115
   272   111
   488   128
   474   159
   414   155
   369   156
   349   146
   279   170
    62   165
   142   169
   212   199
   138   207
    72   200
   152   229
   191   233
   374   227
   457   223
   345   276
   401   316
   298   297
   295   329
   239   268
   245   303
   225   291
   211   264
   115   288
    68   308
    81   280
    92   373
   226   385
   355   466
   466   441
   250   479
   161   469
    93   428
    80   442
    66   454
    91   473
   169   126
   199   122
   375    29
   420    55
   477    27];
   
nn = 40; %20;

if 0
    j=1;
    [res, z, zz, a, b, model] = Hiroshi(fname1,0,[pos(j,2)-nn pos(j,1)-nn; pos(j,2)+nn pos(j,1)+nn]);
    save HiroshiModel model
    eval(['save ' outfile 'M' mint2str(j,3) ' res z zz a b']);
    cnt = 1; bld = 1;
    for jj=1:size(z,3)
        subplot(4,5,cnt+5+5*floor((cnt-1)/5));
        mim(zz(:,:,jj));
        if ismember(jj,res(:,1))
            title({['\theta = ' mint2str(res(res(:,1)==jj,4)) '°']; ...
                ['\phi = ' mint2str(res(res(:,1)==jj,5)) '°']},'fontsize',12,'fontname','times');
        end
        subplot(4,5,cnt+5*floor((cnt-1)/5));
        mim(z(:,:,jj));
        title(num2str(jj),'fontsize',12);
        cnt = cnt+1;
        if cnt==11 | jj==ceil(size(z,3))
            eval(['print -dpng ' outfile 'M' mint2str(j,3) '-' mint2str(bld,3)]);
            bld = bld+1;
            cnt = 1;
        end
    end
    jstart = 2;
else
    load HiroshiModel
    jstart = 7;
end

for j=jstart:size(pos,1)
    [res, z, zz, a, b] = Hiroshi(fname1,0,[pos(j,2)-nn pos(j,1)-nn; pos(j,2)+nn pos(j,1)+nn],model);
    eval(['save ' outfile 'M' mint2str(j,3) ' res z zz a b']);
    cnt = 1; bld = 1; tst = 0;
    for jj=1:size(z,3)
        subplot(4,5,cnt+5+5*floor((cnt-1)/5));
        mim(zz(:,:,jj));
        if ismember(jj,res(:,1))
            title({['\theta = ' mint2str(res(res(:,1)==jj,4)) '°']; ...
                ['\phi = ' mint2str(res(res(:,1)==jj,5)) '°']},'fontsize',12,'fontname','times');
        else
            tst = tst+1;
        end
        subplot(4,5,cnt+5*floor((cnt-1)/5));
        mim(z(:,:,jj));
        title(num2str(jj),'fontsize',12);
        cnt = cnt+1;
        if cnt==11 | jj==ceil(size(z,3))
            if tst<10
                eval(['print -dpng ' outfile 'M' mint2str(j,3) '-' mint2str(bld,3)]);
            end
            bld = bld+1;
            cnt = 1;
            tst = 0;
        end
    end
end


