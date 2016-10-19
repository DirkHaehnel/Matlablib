% readTiff --- reads (stacked) .tif file as unsigned 16Bit integer
% 
% INPUT:
% (optional) Filename (String)
% (optional) First Image of Stack (Integer)
% (optional) Last Image of Stack (Integer)
% (optional) 1 if .tif should be saved as .mat
% 
% OUTPUT:
% FinalImage (Matrix), usage: FinalImage(row,column,stack number)
% FinalImage(1,1,1): Upper left corner of first image
% 
% Jan Thiart 26.9.2012 (from
% http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/)

function [FinalImage] = readTiff(varargin)

if nargin == 0
    [filename,pathname]=uigetfile('*.tif', 'select image file');
    varargin{1}=[pathname,filename];
end    
    
FileTif=varargin{1};
matfilename = regexprep(FileTif, '.tif', '.mat');

if (exist(matfilename,'file') > 0)
    % load existing file
    load(matfilename,'FinalImage');
    
else
    InfoImage=imfinfo(FileTif);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;

        if nargin >= 3
            ImageStart=varargin{2};
            ImageEnd=varargin{3};
            NumberImages=abs(ImageEnd-ImageStart)+1;
        elseif nargin == 2
            ImageStart=varargin{2};
            ImageEnd=ImageStart;
            NumberImages=1;
        elseif nargin == 1
            NumberImages=length(InfoImage);
            ImageStart=1;
            ImageEnd=NumberImages;
        else
            error('Input method: readTiff(Filename,Start of Stack,End of Stack) or \n readTiff(Filename)');
        end

    % read unsigned 16Bit Integer - change if necessary!
    FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
    FileID = tifflib('open',FileTif,'r');
    rps = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);

        for i=ImageStart:ImageEnd
        % for i=1:NumberImages
           tifflib('setDirectory',FileID,i);
           % Go through each strip of data.
           rps = min(rps,mImage);
           for r = 1:rps:mImage
              row_inds = r:min(mImage,r+rps-1);
              stripNum = tifflib('computeStrip',FileID,r);
              FinalImage(row_inds,:,i-ImageStart+1) = tifflib('readEncodedStrip',FileID,stripNum);
              r
           end
           i
        end
        if (nargin == 4 && varargin{4} == 1)
            save(matfilename, 'FinalImage');
        end
        tifflib('close',FileID);
end