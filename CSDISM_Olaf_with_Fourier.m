% spinning disk confocal ISM
% from Joerg, modified by Olaf July 9th 2012

% "mag" gives the magnification of the original image, and thus the
% accuracy with which the centers of the Spin Disk Pinholes can be
% determined
% "seq" is the sequence length, a property of the raw data

clear all
warning off all
close all
dirname = 'Z:\TIRF - Spinning Disc\2012_12_16 test\'; %path where the raw data is

if 0 % conversion
    fnames = dir([dirname '*.tif']);
    for j=1:length(fnames)
        eval(['img = ReadTiff_CP(''' dirname fnames(j).name ''');']);
    end
end

if 0 % fitting the patterns
	%load([dirname 'Reference 100x_Spin Disc 488_0008.mat']);
    eval(['img = ReadTiff_CP(uigetfile);']);%eval(['img = ReadTiff_CP(''' dirname 'ref.tif' ''');']);  %filename of the reference raw data
    img = double(img);
    seq = 250;
    mag = 10;
    
    
    tmp=mean(img(:,:,2:seq:end),3);
    mim(tmp);
    
    [a,b] = PickVec;
    [mx,my,wx,wy] = Gauss2D(tmp(a,b));
    ww = mag*ceil(1.5*mean([wx wy]));
    wv = -ww:1:ww;
    [x,y] = meshgrid(wv,wv);
    mask = exp(-2*((x./mag).^2+(y./mag).^2)/mean([wx,wy])^2);
    
   %deleted the actual pattern fitting since STORM does that now
   
    eval(['save ''' dirname 'CSDISMReference488''' '   seq mask ww wx wy wv mag '])
end

if 1 % image computation
%     tmp = dir([dirname 'reference*.mat']);
%     for j=1:length(tmp)
%         fnames{j} = tmp(j).name(1:end-4);
%     end
%    jj = length(tmp);
    clear fnames;
    jj = 0; cnt = 1;
    tmp = dir([dirname 'Spin Disc 488_0013.tif']);  %filename of the raw data
    for j=1:length(tmp)
        if isempty(strfind(tmp(j).name,'final')) && ~(exist([dirname tmp(j).name(1:end-4) '_final2.mat'],'file')==2) %try to add 'file' to exist.
            fnames{jj+cnt} = tmp(j).name(1:end-4);
            cnt = cnt+1;
        end
    end
    
    load([dirname 'CSDISMReference488']);
    
    load([dirname 'coordinates488']);
    
    load('weight488.mat');
    
    for jf=1:length(fnames)
        %load([dirname '' fnames{jf} '.mat' '']);
        eval(['img = readTiff(''' dirname fnames{jf} '.tif' ''');']);

        if size(img,3)>=seq
            img = double(img);
            ww =floor(size(mask,1)/2);
            wi = (-ww:ww);
         
            res=zeros(2*mag*size(img,1),2*mag*size(img,2));
            %res = zeros(size(xj,1),size(xj,2));
            
            %make image for comparison
            im0=SpreadMatrix(mean(img(:,:,1:floor(size(img,3)/seq)*seq),3),2*mag);
            
            
            
            
            for j = 1:seq %go through the whole sequence
                     
                %magnify the original image
                tmp = SpreadMatrix(mean(img(:,:,j+1:seq:end),3),mag);
                               
                %Define Background
                img_tmp=mean(img(:,:,j+1:seq:end),3);
                [h, bin] = mHist(img_tmp(:),1:max(img_tmp(:)));  %changed tmp to img_tmp
                bck = mean(bin(h==max(h)));
                %c1=im(a,b); c1=c1(:); c1=sort(c1); bck1=mean(c1(1:ceil(end/100)));
                 
                %Remove coordinates which are too close to the boundaries
                ind = xc{j}<=round(ww) | xc{j}>=size(tmp,2)-round(ww) | yc{j}<=round(ww) | yc{j}>=size(tmp,1)-round(ww);
                xc{j}(ind) = []; yc{j}(ind) = [];
                
                %d=Disk(ww);
                
                
                %make higher res image, going through the coordinates of
                %all spinning disk points.
                for k=1:length(xc{j}) %add data from every spin disk pinhole                    
                    res(2*yc{j}(k)+wi,2*xc{j}(k)+wi) = res(2*yc{j}(k)+wi,2*xc{j}(k)+wi) + (tmp(yc{j}(k)+wi,xc{j}(k)+wi)) - bck;                    
                end
                
                %mim(cat(3,im0,res))
            end
            
            
            
            make_small_res;
            res=res_small;
            make_small_im0;
            im0=im0_small;
            
            
            %Fourier stuff
            
            
           
            %Fourier frequency cutoff above 2*k_max and smooth transition
            [x,y]=meshgrid(-511.5:511.5);
            filter = exp(-(x.^2+y.^2)/90^2/2);   %before181
            filter(filter>0.134)=0.134;  %take everything under 181(2*k_max), then smooth transition  , before 0.68
            %filter(filter<0.001)=0.00;
            filter=filter/max(filter(:));
            %cutoff = Disk(442.5);   %Disk(221.5);%Disk(442.5);
            %cutoff_big = zeros(1024,1024);
            %cutoff_big(70:end-69,70:end-69) = cutoff;   %(291:end-290,291:end-290)=d;%

            res_f = abs(ifft2(ifftshift(fftshift(fft2(res)).*weight.*filter)));

            tmp = abs(fftshift(fft2(res_f)));
            qfac = sum(sum((x.^2+y.^2).*tmp))/sum(sum(tmp));

            %mim(cat(3,im0,res,res_f))
            eval(['save ''' dirname fnames{jf} '_final2.mat'' im0 res res_f'])
            save([dirname fnames{jf} '_resf.txt'],'res_f','-ascii');
            save([dirname fnames{jf} '_im0.txt'],'im0','-ascii');
            save([dirname fnames{jf} '_res.txt'],'res','-ascii');
            
        end
    end
end


