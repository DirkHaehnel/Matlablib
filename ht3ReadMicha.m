% HydraHarp File Access Demo in Matlab


% This script reads a binary HydraHarp T3 mode data file (*.ht3)
% and displays the headers. Actual T3 event records are written
% to an output file [filename].out 
% We do not keep them in memory because of the huge amount of memory
% this would take in case of large files. The best is to process
% the data on the fly and keep only the results.


% Works with file format version 1.0 only!
% Tested with Matlab 6.
% Peter Kapusta, PicoQuant GmbH, June 2008
% This is demo code. Use at your own risk. No warranties.
% Make sure you have enough memory when loading large files!


clear all;
clc;

[filename, pathname]=uigetfile('*.ht3', 'T3 mode data:', 0, 0);
fid=fopen([pathname filename]);


% fprintf(1,'\n=========================================================================== \n');
% fprintf(1,'  Content of %s : \n', strcat(pathname, filename));
% fprintf(1,'=========================================================================== \n');
% fprintf(1,'\n');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % ASCII file header
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Ident = char(fread(fid, 16, 'char'));
% fprintf(1,'               Ident: %s\n', Ident);
% 
% 
% FormatVersion = deblank(char(fread(fid, 6, 'char')'));
% fprintf(1,'      Format version: %s\n', FormatVersion);
% 
% 
% if not(strcmp(FormatVersion,'1.0'))
%    fprintf(1,'\n\n      Warning: This program is for version 1.0 only. Aborted.');
%    STOP;
% end;
% 
% 
% CreatorName = char(fread(fid, 18, 'char'));
% fprintf(1,'        Creator name: %s\n', CreatorName);
% 
% 
% CreatorVersion = char(fread(fid, 12, 'char'));
% fprintf(1,'     Creator version: %s\n', CreatorVersion);
% 
% 
% FileTime = char(fread(fid, 18, 'char'));
% fprintf(1,'    Time of creation: %s\n', FileTime);
% 
% 
% CRLF = char(fread(fid, 2, 'char'));
% 
% 
% Comment = char(fread(fid, 256, 'char'));
% fprintf(1,'             Comment: %s\n', Comment);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Binary file header
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % The binary file header information is indentical to that in HHD files.
% % Note that some items are not meaningful in the time tagging modes
% % therefore we do not output them. 
% 
% 
% NumberOfCurves = fread(fid, 1, 'int32');
% 
% 
% BitsPerRecord = fread(fid, 1, 'int32');
% fprintf(1,'       Bits / Record: %d\n', BitsPerRecord);
% 
% 
% ActiveCurve = fread(fid, 1, 'int32');
% 
% 
% MeasurementMode = fread(fid, 1, 'int32');
% fprintf(1,'    Measurement Mode: %d\n', MeasurementMode);
% 
% 
% SubMode = fread(fid, 1, 'int32');
% fprintf(1,'            Sub-Mode: %d\n', SubMode);
% 
% 
% Binning = fread(fid, 1, 'int32');
% fprintf(1,'             Binning: %d\n', Binning);
% 
% 
% Resolution = fread(fid, 1, 'double');
% fprintf(1,'          Resolution: %f ps\n', Resolution);
% 
% 
% Offset = fread(fid, 1, 'int32');
% fprintf(1,'              Offset: %d\n', Offset);
% 
% 
% Tacq = fread(fid, 1, 'int32');
% fprintf(1,'    Acquisition Time: %d ms \n', Tacq);
% 
% 
% StopAt = fread(fid, 1, 'int32');
% StopOnOvfl = fread(fid, 1, 'int32');
% Restart = fread(fid, 1, 'int32');
% DispLinLog = fread(fid, 1, 'int32');
% DispTimeAxisFrom = fread(fid, 1, 'int32');
% DispTimeAxisTo = fread(fid, 1, 'int32');
% DispCountAxisFrom = fread(fid, 1, 'int32');
% DispCountAxisTo = fread(fid, 1, 'int32');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% for i = 1:8
% DispCurveMapTo(i) = fread(fid, 1, 'int32');
% DispCurveShow(i) = fread(fid, 1, 'int32');
% end;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% for i = 1:3
% ParamStart(i) = fread(fid, 1, 'float');
% ParamStep(i) = fread(fid, 1, 'float');
% ParamEnd(i) = fread(fid, 1, 'float');
% end;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% RepeatMode = fread(fid, 1, 'int32');
% RepeatsPerCurve = fread(fid, 1, 'int32');
% RepatTime = fread(fid, 1, 'int32');
% RepeatWaitTime = fread(fid, 1, 'int32');
% ScriptName = char(fread(fid, 20, 'char'));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %          Hardware information header 
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% fprintf(1,'-------------------------------------\n'); 
% 
% 
% HardwareIdent = char(fread(fid, 16, 'char'));
% fprintf(1,' Hardware Identifier: %s\n', HardwareIdent);
% 
% 
% HardwarePartNo = char(fread(fid, 8, 'char'));
% fprintf(1,'Hardware Part Number: %s\n', HardwarePartNo);    
%     
% HardwareSerial = fread(fid, 1, 'int32');
% fprintf(1,'    HW Serial Number: %d\n', HardwareSerial);
% 
% 
% nModulesPresent = fread(fid, 1, 'int32');
% fprintf(1,'     Modules present: %d\n', nModulesPresent);
% 
% 
% for i=1:10
% ModelCode(i) = fread(fid, 1, 'int32');    
% VersionCode(i) = fread(fid, 1, 'int32');
% end;
% for i=1:nModulesPresent
% fprintf(1,'      ModuleInfo[%02d]: %08x %08x\n', i-1, ModelCode(i), VersionCode(i));
% end;
% 
% 
% BaseResolution = fread(fid, 1, 'double');
% fprintf(1,'      BaseResolution: %f\n', BaseResolution);
% 
% 
% InputsEnabled = fread(fid, 1, 'ubit64');
% fprintf(1,'      Inputs Enabled: %x\n', InputsEnabled); %actually a bitfield
% 
% 
% InpChansPresent  = fread(fid, 1, 'int32');
% fprintf(1,' Input Chan. Present: %d\n', InpChansPresent);
% 
% 
% RefClockSource  = fread(fid, 1, 'int32');
% fprintf(1,'      RefClockSource: %d\n', RefClockSource);
% 
% 
% ExtDevices  = fread(fid, 1, 'int32');
% fprintf(1,'    External Devices: %x\n', ExtDevices); %actually a bitfield
% 
% 
% MarkerSettings  = fread(fid, 1, 'int32');
% fprintf(1,'     Marker Settings: %x\n', MarkerSettings); %actually a bitfield
% 
% 
% SyncDivider = fread(fid, 1, 'int32');
% fprintf(1,'        Sync divider: %d \n', SyncDivider);
% 
% 
% SyncCFDLevel = fread(fid, 1, 'int32');
% fprintf(1,'      Sync CFD Level: %d mV\n', SyncCFDLevel);
% 
% 
% SyncCFDZeroCross = fread(fid, 1, 'int32');
% fprintf(1,'  Sync CFD ZeroCross: %d mV\n', SyncCFDZeroCross);
% 
% 
% SyncOffset = fread(fid, 1, 'int32');
% fprintf(1,'         Sync Offset: %d\n', SyncOffset);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %          Channels' information header 
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% for i=1:InpChansPresent
% InputModuleIndex(i) = fread(fid, 1, 'int32');    
% InputCFDLevel(i) = fread(fid, 1, 'int32');
% InputCFDZeroCross(i) = fread(fid, 1, 'int32');    
% InputOffset(i) = fread(fid, 1, 'int32');
% 
% 
% fprintf(1,'\n-------------------------------------\n');
% fprintf(1,'Input Channel No. %d\n', i-1);
% fprintf(1,'-------------------------------------\n');
% fprintf(1,'  Input Module Index: %d\n', InputModuleIndex(i));
% fprintf(1,'     Input CFD Level: %d mV\n', InputCFDLevel(i));
% fprintf(1,' Input CFD ZeroCross: %d mV\n', InputCFDZeroCross(i));
% fprintf(1,'        Input Offset: %d\n', InputOffset(i));
% end;
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %                Time tagging mode specific header
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% fprintf(1,'\n-------------------------------------\n');
% for i=1:InpChansPresent
% InputRate(i) = fread(fid, 1, 'int32');    
% fprintf(1,'     Input Rate [%02d]: %d\n', i-1, InputRate(i));
% end;
% 
% 
% fprintf(1,'\-------------------------------------\n');
% 
% 
% SyncRate = fread(fid, 1, 'int32');
% fprintf(1,'           Sync Rate: %d Hz\n', SyncRate);
% 
% 
% StopAfter = fread(fid, 1, 'int32');
% fprintf(1,'          Stop After: %d ms \n', StopAfter);
% 
% 
% StopReason = fread(fid, 1, 'int32');
% fprintf(1,'         Stop Reason: %d\n', StopReason);
% 
% 
% ImgHdrSize = fread(fid, 1, 'int32');
% fprintf(1,' Imaging Header Size: %d bytes\n', ImgHdrSize);
% 
% 
% nRecords = fread(fid, 1, 'uint64');
% fprintf(1,'   Number of Records: %d\n', nRecords);
% 
% 
% % Special header for imaging. How many of the following ImgHdr array elements
% % are actually present in the file is indicated by ImgHdrSize above. 
% % Storage must be allocated dynamically if ImgHdrSize other than 0 is found.
%  
% ImgHdr = fread(fid, ImgHdrSize, 'int32');  % You have to properly interpret ImgHdr if you want to generate an image 
% 
% 
% 
% % The header section end after ImgHdr. Following in the file are only event records. 
% % How many of them actually are in the file is indicated by nRecords in above.
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %  This reads the T3 mode event records
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
cnt_Ofl=0; cnt_M=0; cnt_Err=0;  % just counters
% syncperiod = 1E9/SyncRate;      % in nanoseconds
% fprintf(1,'\n   Sync Period = %5.4f ns\n',syncperiod);

syncperiod = 12.5;


OverflowCorrection = 0;
T3WRAPAROUND=1024;



outfile = [pathname filename(1:length(filename)-4) '.out'];
fpout = fopen(outfile,'W');
fprintf(1,'\nWriting data to %s', outfile);
fprintf(1,'\nThis may take a while...');



% fprintf(fpout,'\n----------------------------------------------------------------------');
% fprintf(fpout,'\n Record#  T3record  true_nSync  Chan  dtime   true_time/ns');
% fprintf(fpout,'\n----------------------------------------------------------------------');



% for i=1:Records
for i=1:1e6


T3Record = fread(fid, 1, 'ubit32');     % all 32 bits:
  
%   +-------------------------------+  +-------------------------------+ 
%   |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|  |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|
%   +-------------------------------+  +-------------------------------+  
  
nsync = bitand(T3Record,1023);       % the lowest 10 bits:
  
%   +-------------------------------+  +-------------------------------+ 
%   | | | | | | | | | | | | | | | | |  | | | | | | |x|x|x|x|x|x|x|x|x|x|
%   +-------------------------------+  +-------------------------------+  
  
dtime = bitand(bitshift(T3Record,-10),32767);   % the next 15 bits:
%   the dtime unit depends on "Resolution" that can be obtained from header



%   +-------------------------------+  +-------------------------------+ 
%   | | | | | | | |x|x|x|x|x|x|x|x|x|  |x|x|x|x|x|x| | | | | | | | | | |
%   +-------------------------------+  +-------------------------------+


channel = bitand(bitshift(T3Record,-25),63);   % the next 6 bits:


%   +-------------------------------+  +-------------------------------+ 
%   | |x|x|x|x|x|x| | | | | | | | | |  | | | | | | | | | | | | | | | | |
%   +-------------------------------+  +-------------------------------+


special = bitand(bitshift(T3Record,-31),1);   % the last bit:


%   +-------------------------------+  +-------------------------------+ 
%   |x| | | | | | | | | | | | | | | |  | | | | | | | | | | | | | | | | |
%   +-------------------------------+  +-------------------------------+


% fprintf(fpout,'\n%7u   %08x   ',i-1,T3Record); 
  
if special == 0   % this means a regular input channel
   true_nSync = OverflowCorrection + nsync;
   %  one nsync time unit equals to "syncperiod" which can be calculated from "SyncRate"
   truetime = true_nSync * syncperiod;
   
   fprintf(fpout,'%9u  %2u  %6u  %f\n', true_nSync, channel, dtime, truetime);
   
   
else    % this means we have a special record
    
    if channel == 63  % overflow of nsync occured
    % fprintf(fpout,'   OFL');
        OverflowCorrection = OverflowCorrection + T3WRAPAROUND;
    cnt_Ofl=cnt_Ofl+1;
    end;
    
    if (channel>=1)&(channel<=15)  % these are markers
    true_nSync = OverflowCorrection + nsync;
    %  one nsync time unit equals to "syncperiod" which can be calculated from "SyncRate"
    truetime = true_nSync * syncperiod;
    cnt_M=cnt_M+1;
    % fprintf(fpout,'%9u M%2u  %6u  %f', true_nSync, channel, dtime, truetime);
    end;    
    
end;


end;


fprintf(fpout,'\n');
fclose(fpout);
fclose(fid);


fprintf(1,'Ready!  \n\n');
fprintf(1,'\nStatistics obtained from the data:\n');
fprintf(1,'\nLast True Sync = %-14.0f, Last t = %14.3f ns,',true_nSync, truetime);
fprintf(1,'\n%i overflows, %i markers.',cnt_Ofl,cnt_M);
fprintf(1,'\n');