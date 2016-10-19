% spinning disk confocal ISM

%clear all
%close all
dirname = 'Z:\TIRF - Spinning Disc\2012-06-15 Olaf test\';

if 0 % conversion
    fnames = dir([dirname '*.tif']);
    for j=1:length(fnames)
        eval(['img = ReadTiff_CP(''' dirname fnames(j).name ''');']);
    end
end

if 0 % fitting the patterns
	%load([dirname 'Reference 100x_Spin Disc 488_0008.mat']);
    eval(['img = ReadTiff_CP(''' dirname 'Reference 488_Spin Disc 488_0001.tif' ''');']);
    img = double(img);
    seq = 250;
    delta = 0.1;
    max_neigbour_dist = 15/delta;
    
    mim(mean(img(:,:,1:seq:end),3));
    [xi,yi] = meshgrid(1:delta:size(img,1),1:delta:size(img,2));
    [a,b] = PickVec;
    [mx,my,wx,wy] = Gauss2D(img(a,b,1));
    ww = ceil(1.5*mean([wx wy]));
    wv = -ww:delta:ww;
    [x,y] = meshgrid(wv,wv);
    mask = exp(-2*(x.^2+y.^2)/mean([wx,wy])^2);
    [x,y] = meshgrid(1:size(img,1),1:size(img,2));
    for j = 1:seq
        tmp = interp2(x,y,mean(img(:,:,j:seq:end),3),xi,yi);
        [err, bim, cim, sim, xc{j}, yc{j}, ~, ~, ~, ~, imm] = FindPattern(tmp,mask,[],[],[],[],'cim>0.5*max(cim(:))');
        mim(cat(3,tmp,imm));
        eval(['print -dpng -r100 tmp' mint2str(j,4)])
    end
    eval(['save ''' dirname 'CSDISMReference''' ' xc yc seq delta mask ww wx wy wv xi yi x y'])
end

if 1 % image computation
%     tmp = dir([dirname 'reference*.mat']);
%     for j=1:length(tmp)
%         fnames{j} = tmp(j).name(1:end-4);
%     end
%    jj = length(tmp);
    clear fnames;
    jj = 0; cnt = 1;
    tmp = dir([dirname '*0003_Z_*.tif']);
    for j=1:length(tmp)
        if isempty(strfind(tmp(j).name,'final')) && ~(exist([dirname tmp(j).name(1:end-4) '_final.mat'])==2)
            fnames{jj+cnt} = tmp(j).name(1:end-4);
            cnt = cnt+1;
        end
    end
   
    load([dirname 'CSDISMReference']);
    for jf=1:length(fnames)
        %load([dirname '' fnames{jf} '.mat' '']);
        eval(['img = ReadTiff_CP(''' dirname fnames{jf} '.tif' ''');']);
        if size(img,3)>=seq
            img = double(img);
            ww = 8;
            wi = round((-ww:delta:ww)/delta);
            apod = Disk(ww/delta).*exp(-2*(wi'.^2*ones(1,length(wi))+ones(length(wi),1)*wi.^2)/(ww/delta).^2);
            apod = apod/mean(apod(:));
            ws = round((-ww:2*delta:ww)/delta);
            phi = 0:0.01:2*pi;
            res = zeros(size(xi,1),size(xi,2));
            im0 = interp2(x,y,mean(img(:,:,1:floor(size(img,3)/seq)*seq),3),xi,yi);
            
            for j = 1:seq
                tmp = interp2(x,y,mean(img(:,:,j:seq:end),3),xi,yi);
                [h, bin] = mHist(tmp(:),1:max(tmp(:)));
                bck = mean(bin(h==max(h)));
                ind = xc{j}<=round(ww/delta) | xc{j}>=size(tmp,2)-round(ww/delta) | yc{j}<=round(ww/delta) | yc{j}>=size(tmp,1)-round(ww/delta);
                xc{j}(ind) = []; yc{j}(ind) = [];
                for k=1:length(xc{j})
                    res(yc{j}(k)+ws/2,xc{j}(k)+ws/2) = res(yc{j}(k)+ws/2,xc{j}(k)+ws/2) + interp2(wi,wi',tmp(yc{j}(k)+wi,xc{j}(k)+wi)-bck,ws,ws');
                end
                mim(cat(3,im0(3*ww/delta:end-3*ww/delta,3*ww/delta:end-3*ww/delta),res(3*ww/delta:end-3*ww/delta,3*ww/delta:end-3*ww/delta)))
            end
            im0 = im0(3*ww/delta:end-3*ww/delta,3*ww/delta:end-3*ww/delta);
            res = res(3*ww/delta:end-3*ww/delta,3*ww/delta:end-3*ww/delta,:);
            
            eval(['save ''' dirname fnames{jf} '_final.mat'' im0 res'])
            eval(['print -dpng -r300 ''' dirname fnames{jf} '_final.png'''])
        end
    end
end

