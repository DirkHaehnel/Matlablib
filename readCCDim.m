function out = readCCDim(filename,infor,fn,mask);
% This function load fn-th frame of spe file
% filename is the location of the spe file you want to load
% fn is number of frame you want to read
% mask is the data contained in the mat file
% if infor=1 to know data length
% if infor=2 to get the matrix
% type=1:SPE
% type=2:sif
% type=3:tif
% type=4:mat

if strcmpi(filename(end-2:end),'spe')==1
    type=1;
    skp=4100;
elseif strcmpi(filename(end-2:end),'sif')==1
    type=2;
elseif strcmpi(filename(end-2:end),'tif')==1
    type=3;
elseif strcmpi(filename(end-2:end),'mat')==1
    type=4;
    if nargin==3
        load(filename)
    end
end


if type==2
    fid = fopen(filename,'r');
    F = fread(fid,4000);
    s = char(F');
    numb=findstr(s,'65538');
    len=str2num(s(numb(end-1):numb(end)-1));
    ysize=len(1,3);
    xsize=len(1,4);
    datalength=len(1,6);
    s(numb(end):numb(end)+31);
    numb2=findstr(s,num2str(xsize));
    skp=numb2(end)+11+datalength*2;
    s(skp-10:skp+10);

elseif type==1
    spefile = fopen(filename,'r','ieee-le');
    [header,count]=fread(spefile,1100,'uint16');
    header=header';
    fclose(spefile);
    %%%%%%%%%%check header%%%%%%%%%%%%%%%%%%5
    % fid=fopen('header.dat','w');
    % for i=1:4100
    % fprintf(fid, '%d\t%d\r\n', i, header(i));
    % end
    % fclose(fid);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ysize=header(329);
    xsize=header(22);
    datalength=header(724);
    datatype=header(55); %datatype 3*unit16, 2:int16, 1:int32, 0:float32 
    clear header count;
end

%if fn > datalength
%    errordlg('Frame number bigger than data length...')  
%end   

if infor == 1
    if type==3
        datalength=length(imfinfo(filename));%1500;%
    elseif type==4
        datalength=size(mask,3);
    end
    out=datalength;
    imFrame=[];
    
elseif infor == 2
    if type==1
        spefile = fopen(filename,'r','ieee-le');
        if datatype==3;
            st=fseek(spefile, xsize*ysize*(fn-1)*2+skp,-1);
            if st==0;
                imFrame=fread(spefile, xsize*ysize, 'uint16');
            else
                'error!!!'
                return
            end
        elseif datatype ==2;
            st=fseek(spefile, xsize*ysize*(fn-1)*2+skp,-1);
            if st==0;
                imFrame=fread(spefile, xsize*ysize, 'int16');%'*int16'
            else
                'error!!!'
                return
            end
        elseif datatype ==1
            st=fseek(spefile, xsize*ysize*(fn-1)*2+skp,-1);
            if st==0;
                imFrame=fread(spefile, xsize*ysize, 'int32');
            else
                'error!!!'
                return
            end
        elseif datatype==0
            st=fseek(spefile, xsize*ysize*(fn-1)*2+skp,-1);
            if st==0;
                imFrame=fread(spefile, xsize*ysize, 'float32');
            else
                'error!!!'
                return
            end
            clear dam
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fclose(spefile);
        
        %%%%%%remove header and restructure data as 2D matrix%%%%%%%
        imFrame=reshape(imFrame, xsize,ysize);
        imFrame=imFrame';
        
    elseif type==2
        spefile = fopen(filename,'r');   %'ieee-le'
        st=fseek(spefile, skp+(xsize*ysize*(fn-1)*4),-1);
        imFrame=fread(spefile, xsize*ysize,'float32'); 
        imFrame=reshape(imFrame, xsize,ysize);
        imFrame=flipud(imFrame');
        fclose(spefile);
        %imagesc(imFrame'); colormap hot
    elseif type==3
        imFrame=double(imread(filename,fn));
    elseif type==4
        imFrame=mask(:,:,fn);
    end
    out=imFrame;
end
