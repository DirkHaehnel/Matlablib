% spinning disk confocal ISM
% from Joerg, modified by Olaf July 9th 2012

% "mag" gives the magnification of the original image, and thus the
% accuracy with which the centers of the Spin Disk Pinholes can be
% determined
% "seq" is the sequence length, a property of the raw data

%clear all
warning off all
%close all
dirname = 'Z:\TIRF - Spinning Disc\2012-8-5 Olaf blue red\'; %path where the raw data is

if 0 % conversion
    fnames = dir([dirname '*.tif']);
    for j=1:length(fnames)
        eval(['img = ReadTiff_CP(''' dirname fnames(j).name ''');']);
    end
end

if 0 % fitting the patterns
	%load([dirname 'Reference 100x_Spin Disc 488_0008.mat']);
    eval(['img = ReadTiff_CP(uigetfile);']);  %filename of the reference raw data
    img = double(img);
    seq = 500;
    mag = 10;
    
    mim(mean(img(:,:,1:seq:end),3));
    
    [a,b] = PickVec;
    [mx,my,wx,wy] = Gauss2D(img(a,b,1));
    ww = mag*ceil(1.5*mean([wx wy]));
    wv = -ww:1:ww;
    [x,y] = meshgrid(wv,wv);
    mask = exp(-2*((x./mag).^2+(y./mag).^2)/mean([wx,wy])^2);
    
    
    tmp=zeros(mag*size(img,1),mag*size(img,2));
    
    for j = 1:1      
    tmp=SpreadMatrix(mean(img(:,:,j+1:seq:end),3),mag); 

        [~, ~, ~, ~, xc1{j}, yc1{j}, ~, ~, ~, ~, imm] = FindPattern_olaf(tmp,mask,[],[],[],[],'cim>0.5*max(cim(:))');

        mim(cat(3,tmp,imm)) %Display the original and the fitted image
        
        %eval(['print -dpng -r100 tmp' mint2str(j,4)])  %Print image to
        %file (for later reference)
        
    end
    eval(['save ''' dirname 'CSDISMReference''' '   seq mask ww wx wy wv mag '])
end

if 1 % image computation
%     tmp = dir([dirname 'reference*.mat']);
%     for j=1:length(tmp)
%         fnames{j} = tmp(j).name(1:end-4);
%     end
%    jj = length(tmp);
    clear fnames;
    jj = 0; cnt = 1;
    tmp = dir([dirname 'cell*647_*_Z*.tif']);  %filename of the raw data
    for j=1:length(tmp)
        if isempty(strfind(tmp(j).name,'final')) && ~(exist([dirname tmp(j).name(1:end-4) '_final.mat'],'file')==2) %try to add 'file' to exist.
            fnames{jj+cnt} = tmp(j).name(1:end-4);
            cnt = cnt+1;
        end
    end
    
    load([dirname 'CSDISMReference']);
    seq=500;
    
    for jf=1:length(fnames)
        %load([dirname '' fnames{jf} '.mat' '']);
        eval(['img = ReadTiff_CP(''' dirname fnames{jf} '.tif' ''');']);

        if size(img,3)>=seq
            img = double(img);
            ww =floor(size(mask,1)/2);
            wi = (-ww:ww);
         
            res=zeros(2*mag*size(img,1),2*mag*size(img,2));
            %res = zeros(size(xj,1),size(xj,2));
            
            %make image for comparison
            im0=SpreadMatrix(mean(img(:,:,1:floor(size(img,3)/seq)*seq),3),2*mag);
           
fnames{jf}
            for j = 1:seq %go through the whole sequence
                     j
                %magnify the original image
                tmp = SpreadMatrix(mean(img(:,:,j+1:seq:end),3),mag);
                               
                %Define Background
                [h, bin] = mHist(tmp(:),1:max(tmp(:)));
                bck = mean(bin(h==max(h)));
                 
                %Remove coordinates which are too close to the boundaries
                ind = xc{j}<=round(ww) | xc{j}>=size(tmp,2)-round(ww) | yc{j}<=round(ww) | yc{j}>=size(tmp,1)-round(ww);
                xc{j}(ind) = []; yc{j}(ind) = [];
                
                %make higher res image, going through the coordinates of
                %all spinning disk points.
                for k=1:length(xc{j}) %add data from every spin disk pinhole                    
                    res(2*yc{j}(k)+wi,2*xc{j}(k)+wi) = res(2*yc{j}(k)+wi,2*xc{j}(k)+wi) + tmp(yc{j}(k)+wi,xc{j}(k)+wi)-bck;                    
                end
                
                mim(cat(3,im0,res))
            end
            
            
            eval(['save ''' dirname fnames{jf} '_final_1.mat'' im0 res'])
            %eval(['print -dpng -r300 ''' dirname fnames{jf} '_final.png'''])
            make_small_res;
            make_small_im0;
            save([dirname fnames{jf} '_res_1.txt'],'res_small','-ascii');
            save([dirname fnames{jf} '_im0_1.txt'],'im0_small','-ascii');
        end
    end
end


