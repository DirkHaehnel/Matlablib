function [y, tcspc, chan, markers, num, overcount, head] = pt3v2read(name, cnts)

% [y, tcspc, chan, markers, num, overcount, head] = pt3v2read(name, cnts)
% (c) Peter Kapusta, Michael Wahl, PicoQuant GmbH 2007, updated May 2007

fid = fopen(name);

%
% The following represents the readable ASCII file header portion 
%

head.Ident = char(fread(fid, 16, 'char')');
head.FormatVersion = deblank(char(fread(fid, 6, 'char')'));

if not(strcmp(head.FormatVersion,'2.0'))
   fprintf(1,'\n\n      Warning: This program is for version 2.0 only. Aborted.');
   fclose(fid);
   return;
end;

head.CreatorName = char(fread(fid, 18, 'char')');
head.CreatorVersion = char(fread(fid, 12, 'char')');
head.FileTime = char(fread(fid, 18, 'char')');
head.CRLF = char(fread(fid, 2, 'char')');
head.CommentField = char(fread(fid, 256, 'char')');

%
% The following is binary file header information
%

head.Curves = fread(fid, 1, 'int32');
head.BitsPerRecord = fread(fid, 1, 'int32');
head.RoutingChannels = fread(fid, 1, 'int32');
head.NumberOfBoards = fread(fid, 1, 'int32');
head.ActiveCurve = fread(fid, 1, 'int32');
head.MeasurementMode = fread(fid, 1, 'int32');
head.SubMode = fread(fid, 1, 'int32');
head.RangeNo = fread(fid, 1, 'int32');
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
head.ScriptName = fread(fid, 20, 'char');

%
% The next is a board specific header
%

head.HardwareIdent = fread(fid, 16, 'char');
head.HardwareVersion = fread(fid, 8, 'char');
head.HardwareSerial = fread(fid, 1, 'int32');
head.SyncDivider = fread(fid, 1, 'int32');
head.CFDZeroCross0 = fread(fid, 1, 'int32');
head.CFDLevel0 = fread(fid, 1, 'int32');
head.CFDZeroCross1 = fread(fid, 1, 'int32');
head.CFDLevel1 = fread(fid, 1, 'int32');
head.Resolution = fread(fid, 1, 'float');

% below is new in format version 2.0

head.RouterModelCode      = fread(fid, 1, 'int32');
head.RouterEnabled        = fread(fid, 1, 'int32');

% Router Ch1
head.RtChan1_InputType    = fread(fid, 1, 'int32');
head.RtChan1_InputLevel   = fread(fid, 1, 'int32');
head.RtChan1_InputEdge    = fread(fid, 1, 'int32');
head.RtChan1_CFDPresent   = fread(fid, 1, 'int32');
head.RtChan1_CFDLevel     = fread(fid, 1, 'int32');
head.RtChan1_CFDZeroCross = fread(fid, 1, 'int32');
% Router Ch2
head.RtChan2_InputType    = fread(fid, 1, 'int32');
head.RtChan2_InputLevel   = fread(fid, 1, 'int32');
head.RtChan2_InputEdge    = fread(fid, 1, 'int32');
head.RtChan2_CFDPresent   = fread(fid, 1, 'int32');
head.RtChan2_CFDLevel     = fread(fid, 1, 'int32');
head.RtChan2_CFDZeroCross = fread(fid, 1, 'int32');
% Router Ch3
head.RtChan3_InputType    = fread(fid, 1, 'int32');
head.RtChan3_InputLevel   = fread(fid, 1, 'int32');
head.RtChan3_InputEdge    = fread(fid, 1, 'int32');
head.RtChan3_CFDPresent   = fread(fid, 1, 'int32');
head.RtChan3_CFDLevel     = fread(fid, 1, 'int32');
head.RtChan3_CFDZeroCross = fread(fid, 1, 'int32');
% Router Ch4
head.RtChan4_InputType    = fread(fid, 1, 'int32');
head.RtChan4_InputLevel   = fread(fid, 1, 'int32');
head.RtChan4_InputEdge    = fread(fid, 1, 'int32');
head.RtChan4_CFDPresent   = fread(fid, 1, 'int32');
head.RtChan4_CFDLevel     = fread(fid, 1, 'int32');
head.RtChan4_CFDZeroCross = fread(fid, 1, 'int32');

% Router settings are meaningful only for an existing router:

%
% The next is a T3 mode specific header
%

head.ExtDevices = fread(fid, 1, 'int32');
head.Reserved1 = fread(fid, 1, 'int32');
head.Reserved2 = fread(fid, 1, 'int32');
head.CntRate0 = fread(fid, 1, 'int32');
head.CntRate1 = fread(fid, 1, 'int32');
head.StopAfter = fread(fid, 1, 'int32');
head.StopReason = fread(fid, 1, 'int32');
head.Records = fread(fid, 1, 'uint32');
head.ImgHdrSize = fread(fid, 1, 'int32');

%Special header for imaging 
if head.ImgHdrSize==15
    head.ImgHdr.Dimensions = fread(fid, 1, 'int32');
    head.ImgHdr.Ident = fread(fid, 1, 'int32');    
    head.ImgHdr.TimePerPixel = fread(fid, 1, 'int32');    
    head.ImgHdr.Acceleration = fread(fid, 1, 'int32');   
    head.ImgHdr.Pattern = fread(fid, 1, 'int32');    
    head.ImgHdr.Reserved = fread(fid, 1, 'int32');    
    head.ImgHdr.X0 = fread(fid, 1, 'float');    
    head.ImgHdr.Y0 = fread(fid, 1, 'float');        
    head.ImgHdr.PixX = fread(fid, 1, 'int32');        
    head.ImgHdr.PixY = fread(fid, 1, 'int32');        
    head.ImgHdr.PixResol = fread(fid, 1, 'float');    
    head.ImgHdr.TStartTo = fread(fid, 1, 'float');        
    head.ImgHdr.TStopTo = fread(fid, 1, 'float');    
    head.ImgHdr.TStartFro = fread(fid, 1, 'float');        
    head.ImgHdr.TStopFro = fread(fid, 1, 'float');    
else
    head.ImgHdr = fread(fid, head.ImgHdrSize, 'uint32');
end

%
%  This reads the T3 mode event records
%

WRAPAROUND=65536;

head.Sync = 1E9/head.CntRate0;   % in nanoseconds

if nargin>1
    fseek(fid, 4*(cnts(1)-1), 0);
    [T3Record num] = fread(fid, cnts(2), 'ubit32'); % all 32 bits:
    %   +-------------------------------+  +-------------------------------+
    %   |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|  |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|
    %   +-------------------------------+  +-------------------------------+

    y = bitand(T3Record,65535);       % the lowest 16 bits:
    %   +-------------------------------+  +-------------------------------+
    %   | | | | | | | | | | | | | | | | |  |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|
    %   +-------------------------------+  +-------------------------------+

    chan = bitand(bitshift(T3Record,-28),15);   % the upper 4 bits:
    %   +-------------------------------+  +-------------------------------+
    %   |x|x|x|x| | | | | | | | | | | | |  | | | | | | | | | | | | | | | | |
    %   +-------------------------------+  +-------------------------------+

    tcspc = bitand(bitshift(T3Record,-16),4095);
    % if chan = 1,2,3 or 4, then these  bits contain the dtime:
    %   +-------------------------------+  +-------------------------------+
    %   | | | | |x|x|x|x|x|x|x|x|x|x|x|x|  | | | | | | | | | | | | | | | | |
    %   +-------------------------------+  +-------------------------------+

    markers = bitand(bitshift(T3Record,-16),15); % where these four bits are markers:
    %   +-------------------------------+  +-------------------------------+
    %   | | | | | | | | | | | | |x|x|x|x|  | | | | | | | | | | | | | | | | |
    %   +-------------------------------+  +-------------------------------+

    y = y + WRAPAROUND*cumsum(markers==0 & chan==15);
    overcount = WRAPAROUND*sum(markers==0 & chan==15);
    
    valid = (chan>4)&(chan~=15)|(chan==15)&(tcspc==0);
    y(valid)       = [];
    tcspc(valid)   = [];
    chan(valid)    = [];
    markers(valid) = [];
else
    y = head;
end

fclose(fid);