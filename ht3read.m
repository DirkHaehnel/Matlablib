function [sync, tcspc, chan, special, num, overcount, head] = ht3read(name, cnts)

% [y, chan, special, num, overcount, head] = ht2read(name, cnts)
% (c) Joerg Enderlein, Michael Wahl, (2008)

fid = fopen(name);

%
% The following represents the readable ASCII file header portion 
%

% head.Ident = char(fread(fid, 16, 'char')');
% head.FormatVersion = deblank(char(fread(fid, 6, 'char')'));
% 
% if not(strcmp(head.FormatVersion,'2.0'))
%    fprintf(1,'\n\n      Warning: This program is for version 2.0 only. Aborted.');
%    fclose(fid);
%    return;
% end;
% 
% head.CreatorName = char(fread(fid, 18, 'char')');
% head.CreatorVersion = char(fread(fid, 12, 'char')');
% head.FileTime = char(fread(fid, 18, 'char')');
% head.CRLF = char(fread(fid, 2, 'char')');
% head.CommentField = char(fread(fid, 256, 'char')');
% 
% %
% % The following is binary file header information
% %
% 
% head.Curves = fread(fid, 1, 'int32');
% head.BitsPerRecord = fread(fid, 1, 'int32');
% head.RoutingChannels = fread(fid, 1, 'int32');
% head.NumberOfBoards = fread(fid, 1, 'int32');
% head.ActiveCurve = fread(fid, 1, 'int32');
% head.MeasurementMode = fread(fid, 1, 'int32');
% head.SubMode = fread(fid, 1, 'int32');
% head.RangeNo = fread(fid, 1, 'int32');
% head.Offset = fread(fid, 1, 'int32');
% head.AcquisitionTime = fread(fid, 1, 'int32');
% head.StopAt = fread(fid, 1, 'int32');
% head.StopOnOvfl = fread(fid, 1, 'int32');
% head.Restart = fread(fid, 1, 'int32');
% head.DispLinLog = fread(fid, 1, 'int32');
% head.DispTimeFrom = fread(fid, 1, 'int32');
% head.DispTimeTo = fread(fid, 1, 'int32');
% head.DispCountFrom = fread(fid, 1, 'int32');
% head.DispCountTo = fread(fid, 1, 'int32');
% 
% for i = 1:8
%     head.DispCurveMapTo(i) = fread(fid, 1, 'int32');
%     head.DispCurveShow(i) = fread(fid, 1, 'int32');
% end;
% 
% for i = 1:3
%     head.ParamStart(i) = fread(fid, 1, 'float');
%     head.ParamStep(i) = fread(fid, 1, 'float');
%     head.ParamEnd(i) = fread(fid, 1, 'float');
% end;
% 
% head.RepeatMode = fread(fid, 1, 'int32');
% head.RepeatsPerCurve = fread(fid, 1, 'int32');
% head.RepeatTime = fread(fid, 1, 'int32');
% head.RepeatWait = fread(fid, 1, 'int32');
% head.ScriptName = fread(fid, 20, 'char');
% 
% %
% % The next is a board specific header
% %
% 
% head.HardwareIdent = fread(fid, 16, 'char');
% head.HardwareVersion = fread(fid, 8, 'char');
% head.HardwareSerial = fread(fid, 1, 'int32');
% head.SyncDivider = fread(fid, 1, 'int32');
% head.CFDZeroCross0 = fread(fid, 1, 'int32');
% head.CFDLevel0 = fread(fid, 1, 'int32');
% head.CFDZeroCross1 = fread(fid, 1, 'int32');
% head.CFDLevel1 = fread(fid, 1, 'int32');
% head.Resolution = fread(fid, 1, 'float');
% 
% % below is new in format version 2.0
% 
% head.RouterModelCode      = fread(fid, 1, 'int32');
% head.RouterEnabled        = fread(fid, 1, 'int32');
% 
% % Router Ch1
% head.RtChan1_InputType    = fread(fid, 1, 'int32');
% head.RtChan1_InputLevel   = fread(fid, 1, 'int32');
% head.RtChan1_InputEdge    = fread(fid, 1, 'int32');
% head.RtChan1_CFDPresent   = fread(fid, 1, 'int32');
% head.RtChan1_CFDLevel     = fread(fid, 1, 'int32');
% head.RtChan1_CFDZeroCross = fread(fid, 1, 'int32');
% % Router Ch2
% head.RtChan2_InputType    = fread(fid, 1, 'int32');
% head.RtChan2_InputLevel   = fread(fid, 1, 'int32');
% head.RtChan2_InputEdge    = fread(fid, 1, 'int32');
% head.RtChan2_CFDPresent   = fread(fid, 1, 'int32');
% head.RtChan2_CFDLevel     = fread(fid, 1, 'int32');
% head.RtChan2_CFDZeroCross = fread(fid, 1, 'int32');
% % Router Ch3
% head.RtChan3_InputType    = fread(fid, 1, 'int32');
% head.RtChan3_InputLevel   = fread(fid, 1, 'int32');
% head.RtChan3_InputEdge    = fread(fid, 1, 'int32');
% head.RtChan3_CFDPresent   = fread(fid, 1, 'int32');
% head.RtChan3_CFDLevel     = fread(fid, 1, 'int32');
% head.RtChan3_CFDZeroCross = fread(fid, 1, 'int32');
% % Router Ch4
% head.RtChan4_InputType    = fread(fid, 1, 'int32');
% head.RtChan4_InputLevel   = fread(fid, 1, 'int32');
% head.RtChan4_InputEdge    = fread(fid, 1, 'int32');
% head.RtChan4_CFDPresent   = fread(fid, 1, 'int32');
% head.RtChan4_CFDLevel     = fread(fid, 1, 'int32');
% head.RtChan4_CFDZeroCross = fread(fid, 1, 'int32');
% 
% % Router settings are meaningful only for an existing router:
% 
% %
% % The next is a T3 mode specific header
% %
% 
% head.ExtDevices = fread(fid, 1, 'int32');
% head.Reserved1 = fread(fid, 1, 'int32');
% head.Reserved2 = fread(fid, 1, 'int32');
% head.CntRate0 = fread(fid, 1, 'int32');
% head.CntRate1 = fread(fid, 1, 'int32');
% head.StopAfter = fread(fid, 1, 'int32');
% head.StopReason = fread(fid, 1, 'int32');
% head.Records = fread(fid, 1, 'uint32');
% head.ImgHdrSize = fread(fid, 1, 'int32');
% 
% %Special header for imaging 
% head.ImgHdr = fread(fid, head.ImgHdrSize, 'int32');

%
%  This reads the HydraHarp binary data
%

WRAPAROUND=1024;

head.Resolution = 2e-3;   % in nanoseconds

if nargin>1
    fseek(fid, 4*(cnts(1)-1), 0);
    [T3Record num] = fread(fid, cnts(2), 'ubit32'); % all 32 bits:
    %   +-------------------------------+  +-------------------------------+
    %   |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|  |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|
    %   +-------------------------------+  +-------------------------------+

    special = bitand(bitshift(T3Record,-31),1);       
    %   +-------------------------------+  +-------------------------------+
    %   |x| | | | | | | | | | | | | | | |  | | | | | | | | | | | | | | | | |
    %   +-------------------------------+  +-------------------------------+

    chan = bitand(bitshift(T3Record,-25),63);   
    %   +-------------------------------+  +-------------------------------+
    %   | |x|x|x|x|x|x| | | | | | | | | |  | | | | | | | | | | | | | | | | |
    %   +-------------------------------+  +-------------------------------+
   
    tcspc = bitand(bitshift(T3Record,-10),32767); 
    %   +-------------------------------+  +-------------------------------+
    %   | | | | | | | |x|x|x|x|x|x|x|x|x|  |x|x|x|x|x|x| | | | | | | | | | |
    %   +-------------------------------+  +-------------------------------+

    sync = bitand(T3Record,1023); 
    %   +-------------------------------+  +-------------------------------+
    %   | | | | | | | | | | | | | | | | |  | | | | | | |x|x|x|x|x|x|x|x|x|x|
    %   +-------------------------------+  +-------------------------------+

    ind = special==1 & chan==63;
    sync = sync + WRAPAROUND*cumsum(ind);
    overcount = WRAPAROUND*sum(ind);
    sync(ind) = [];
    tcspc(ind) = [];
    special(ind) = [];
    chan(ind) = [];
else
    sync = head;
end

fclose(fid);

