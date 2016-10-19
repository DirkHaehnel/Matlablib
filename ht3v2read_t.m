function [sync, tcspc, chan, special, num, head] = ht3v2read_t(name, cnts)
%
%  function [sync, tcspc, chan, special, num, head] = ht3v2read_t(name, cnts)
%
%  This function reads single-photon data from the file 'name'
%
%  If the argument parameter 'cnts' is missing, or has the value 0 it returns just the 
%  header of the file in the output variable 'sync'.
%  
%  If 'cnts' contains two numbers [cnts(1) cnts(2)], the routine proceeds
%  to the position cnts(1) and reads the data until the number of sync
%  events exceeds cnts(2). By this one reads a fixed time period of the
%  photon stream (t = cnts(2) / repetition rate).
%
%  If 'cnts' contains three numbers [cnts(1) cnts(2) cnts(3)], the routine
%  does not read the header, but jumps over the first cnts(3)-bytes that
%  represent the length of the file-header. After that it procceds to the
%  data position cnts(1) before reading the necessary records of data.
%
%  The output variables contain the followig data:
%  sync    : number of the sync events that preceeded this detection event
%  tcspc   : number of the tcspc-bin of the event
%  chan    : number of the input channel of the event (detector-number)
%  special : indicator of the event-type (0: photon; else : virtual photon)
%  num     : counter of the records that were actually read 
%  head    : the file header 
%

if (nargin<2)||isempty(cnts)
    cnts = [0 0 0];
end;

if numel(cnts)<3
    cnts = [cnts 0];
end;

if numel(cnts)<2
    cnts = [0 cnts 0];
end;

if isempty(name)
    fprintf(1,'\n\n      You have to specify a valid file-name. Aborted.\n');
    return;
else    
    fid = fopen(name);
    if fid<1
        fprintf(1,'\n\n      Could not open <%s>. Aborted.\n', name);
        return;
    end
end

if (nargout==6)||(cnts(3)==0)
    
    head = [];
    
    tmp = deblank(char(fread(fid, 16, 'char')'));

    if strcmp(tmp, 'HydraHarp')
        
        head.Ident = tmp;
        head.FormatVersion = deblank(char(fread(fid, 6, 'char')'));

        if not(strcmp(head.FormatVersion,'1.0'))
            fprintf(1,'\n\n      Warning: This program is for version 1.0 only. Aborted.\n');
            fclose(fid);
            return;
        end;

        head.CreatorName = deblank(char(fread(fid, 18, 'char')'));
        head.CreatorVersion = deblank(char(fread(fid, 12, 'char')'));
        head.FileTime = deblank(char(fread(fid, 18, 'char')'));
        fread(fid, 2, 'char');
        head.CommentField = deblank(char(fread(fid, 256, 'char')'));

        %
        % The following is binary file header information
        %

        head.Curves = fread(fid, 1, 'int32');
        head.BitsPerRecord = fread(fid, 1, 'int32');
        head.ActiveCurve = fread(fid, 1, 'int32');
        head.MeasurementMode = fread(fid, 1, 'int32');
        head.SubMode = fread(fid, 1, 'int32');
        head.Binning = fread(fid, 1, 'int32');
        head.Resolution = fread(fid, 1, 'double')/1000;
        head.Offset = fread(fid, 1, 'int32');
        head.AcquisitionTime = fread(fid, 1, 'int32');
        head.StopAt = fread(fid, 1, 'int32');
        head.StopOnOvfl = fread(fid, 1, 'int32');
        head.Restart = fread(fid, 1, 'int32');
        head.DispLinLog = fread(fid, 1, 'int32');
        head.DispTimeFrom = fread(fid, 1, 'int32');
        head.DispTimeTo = fread(fid, 1, 'int32');
        head.DispCountFrom = fread(fid, 1, 'int32');
        head.DispCountTo = fread(fid, 1, 'int32');

        for i = 1:8
            head.DispCurveMapTo(i) = fread(fid, 1, 'int32');
            head.DispCurveShow(i) = fread(fid, 1, 'int32');
        end;

        for i = 1:3
            head.ParamStart(i) = fread(fid, 1, 'float');
            head.ParamStep(i) = fread(fid, 1, 'float');
            head.ParamEnd(i) = fread(fid, 1, 'float');
        end;

        head.RepeatMode = fread(fid, 1, 'int32');
        head.RepeatsPerCurve = fread(fid, 1, 'int32');
        head.RepeatTime = fread(fid, 1, 'int32');
        head.RepeatWait = fread(fid, 1, 'int32');
        head.ScriptName = deblank(char(fread(fid, 20, 'char')'));

        %
        % The next is a board specific header
        %

        head.HardwareIdent = deblank(char(fread(fid, 16, 'char')'));
        head.HardwareVersion = deblank(char(fread(fid, 8, 'char')'));
        head.HardwareSerial = fread(fid, 1, 'int32');
        head.NumModules  = fread(fid, 1, 'int32');

        for i = 1:10
            head.Module(i).Serial = fread(fid, 1, 'int32');
            head.Module(i).Version = fread(fid, 1, 'int32');
        end;

        head.BaseResolution = fread(fid, 1, 'double');
        head.InputsEnabled = fread(fid,1,'int64');
        head.InputsPresent = fread(fid,1,'int32');
        head.RefClock = fread(fid,1,'int32');
        head.ExtDevices = fread(fid,1,'int32');
        head.MarkerSettings = fread(fid,1,'int32');

        head.SyncDivider = fread(fid, 1, 'int32');
        head.SyncCFDLevel = fread(fid, 1, 'int32');
        head.SyncCFDZeroCross = fread(fid, 1, 'int32');
        head.SyncOffset = fread(fid, 1, 'int32');

        for i = 1:head.InputsPresent
            head.InpChannel(i).ModuleIdx = fread(fid, 1, 'int32');
            head.InpChannel(i).CFDLevel = fread(fid, 1, 'int32');
            head.InpChannel(i).CFDZeroCross = fread(fid, 1, 'int32');
            head.InpChannel(i).Offset = fread(fid, 1, 'int32');
        end;

        for i = 1:head.InputsPresent
            head.InpChannel(i).CntRate = fread(fid, 1, 'int32');
        end;

        %
        % The next is a T3 mode specific header
        %

        head.SyncRate = fread(fid, 1, 'int32');
        head.StopAfter = fread(fid, 1, 'int32');
        head.StopReason = fread(fid, 1, 'int32');
        head.ImgHdrSize = fread(fid, 1, 'int32');
        head.Records = fread(fid, 1, 'int64');

        %Special header for imaging

        if head.ImgHdrSize > 0
            head.ImgHdr.Dimensions = fread(fid, 1, 'int32');
            head.ImgHdr.Ident = fread(fid, 1, 'int32');
            if head.ImgHdr.Ident == 1
                head.ImgHdr.PixelTime =  5e6*fread(fid, 1, 'int32');
                head.ImgHdr.Acceleration =  fread(fid, 1, 'int32');
                head.ImgHdr.Pattern =  fread(fid, 1, 'int32');
                fread(fid, 1, 'int32');
                head.ImgHdr.X0 = fread(fid, 1, 'float');
                head.ImgHdr.Y0 = fread(fid, 1, 'float');
                head.ImgHdr.PixX =  fread(fid, 1, 'int32');
                head.ImgHdr.PixY =  fread(fid, 1, 'int32');
                head.ImgHdr.PixelSize =  fread(fid, 1, 'float');
                head.ImgHdr.TStartTo  =  fread(fid, 1, 'float');
                head.ImgHdr.TStopTo   =  fread(fid, 1, 'float');
                head.ImgHdr.TStartFro =  fread(fid, 1, 'float');
                head.ImgHdr.TStopFro  =  fread(fid, 1, 'float');

            else
                head.ImgHdr.data = fread(fid, head.ImgHdrSize-2, 'int32');
            end
        end;
        head.length = ftell(fid);
    else
        head.Resolution = 2e-3;   % in nanoseconds
        head.length = 0;
    end
end


if cnts(2)>0

    if cnts(3)>0
        fseek(fid, (cnts(3)), 'bof');
    end
    if cnts(1)>0
        fseek(fid, 4*(cnts(1)-1), 'cof');
    end
    
    WRAPAROUND=1024;
    
    special = [];
    sync    = 0;
    tcspc   = [];
    chan    = [];

    while (sync(end)<cnts(2))&&(~feof(fid))
        
        [T3Record] = fread(fid, 5E5, 'ubit32');
                
        n_special = bitand(bitshift(T3Record,-31),    1);
        n_chan    = bitand(bitshift(T3Record,-25),   63);
        n_tcspc   = bitand(bitshift(T3Record,-10),32767);        
        n_sync    = bitand(         T3Record     , 1023);
        
        ind     = (n_special==1 & n_chan==63);
        cs_ovl  = sync(end) + WRAPAROUND*cumsum(ind);
        
        special = [special; n_special];       %#ok<AGROW>
        chan    = [chan;    n_chan];          %#ok<AGROW>
        tcspc   = [tcspc;   n_tcspc];         %#ok<AGROW>
        sync    = [sync;    n_sync + cs_ovl]; %#ok<AGROW>

    end;
    
    sync(1) = [];
    ind     = (sync>cnts(2));
        
    if (sum(ind)>0)
            sync(ind)     = [];
            tcspc(ind)    = [];
            special(ind)  = [];
            chan(ind)     = [];
    end;
    
    num = numel(sync);

    if  num>0
        ind    = (special==1 & chan==63);
        sync(ind)     = [];
        tcspc(ind)    = [];
        special(ind)  = [];
        chan(ind)     = [];
    end

else
    sync = head;
end

fclose(fid);
    

