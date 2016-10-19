% HiroshiMaster

clear all

fname1 = {'D:\Daten\Johan\G0_PMA_001\G0_PMA_001.TIF'};
outfile = 'D:\Daten\Johan\G0_PMA_001\G0_PMA_001';

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

pos = [...
    27    90
    36    24
    85    36
   120   107
    58   160
   177    47];

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


