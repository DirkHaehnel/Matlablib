% Demo accessing PicoHarp data files (*.phd) from MATLAB
% PicoHarp 300, Software version 1.0
% Tested with Matlab 5 and 6.
% Peter Kapusta, PicoQuant GmbH, January 2005
% This is demo code. Use at your own risk. No warranties.
% Make sure you have enough memory when loading large files!

clear all;
clc;

[filename, pathname]=uigetfile('*.phd', 'Interactive mode data:', 0, 0);
fid=fopen([pathname filename]);

fprintf(1,'================================================================= \n');
fprintf(1,'  Content of %s : \n', strcat(pathname, filename));
fprintf(1,'================================================================= \n');
fprintf(1,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ASCII file header
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ident = setstr(fread(fid, 16, 'char'));
fprintf(1,'               Ident: %s\n', Ident);

FormatVersion = setstr(fread(fid, 6, 'char'));
fprintf(1,'      Format version: %s\n', FormatVersion);

CreatorName = setstr(fread(fid, 18, 'char'));
fprintf(1,'        Creator name: %s\n', CreatorName);

CreatorVersion = setstr(fread(fid, 12, 'char'));
fprintf(1,'     Creator version: %s\n', CreatorVersion);

FileTime = setstr(fread(fid, 18, 'char'));
fprintf(1,'    Time of creation: %s\n', FileTime);

CRLF = setstr(fread(fid, 2, 'char'));

Comment = setstr(fread(fid, 256, 'char'));
fprintf(1,'             Comment: %s\n', Comment);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Binary file header
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NumberOfCurves = fread(fid, 1, 'int32');
fprintf(1,'    Number of Curves: %d\n', NumberOfCurves);

BitsPerHistoBin = fread(fid, 1, 'int32');
fprintf(1,'     Bits / HistoBin: %d\n', BitsPerHistoBin);

RoutingChannels = fread(fid, 1, 'int32');
fprintf(1,'    Routing channels: %d\n', RoutingChannels);

NumberOfBoards = fread(fid, 1, 'int32');
fprintf(1,'    Number of boards: %d\n', NumberOfBoards);

ActiveCurve = fread(fid, 1, 'int32');
fprintf(1,'        Active Curve: %d\n', ActiveCurve);

MeasurementMode = fread(fid, 1, 'int32');
fprintf(1,'    Measurement Mode: %d\n', MeasurementMode);

SubMode = fread(fid, 1, 'int32');
fprintf(1,'            Sub-Mode: %d\n', SubMode);

RangeNo = fread(fid, 1, 'int32');
fprintf(1,'            Range No: %d\n', RangeNo);

Offset = fread(fid, 1, 'int32');
fprintf(1,'              Offset: %d\n', Offset);

Tacq = fread(fid, 1, 'int32');
fprintf(1,'    Acquisition time: %d ms \n', Tacq);

StopAt = fread(fid, 1, 'int32');
fprintf(1,'             Stop at: %d counts \n', StopAt);

StopOnOvfl = fread(fid, 1, 'int32');
fprintf(1,'    Stop on Overflow: %d\n', StopOnOvfl);

Restart = fread(fid, 1, 'int32');
fprintf(1,'             Restart: %d\n', Restart);

DispLinLog = fread(fid, 1, 'int32');
fprintf(1,'     Display Lin/Log: %d\n', DispLinLog);

DispTimeAxisFrom = fread(fid, 1, 'int32');
fprintf(1,'      Time Axis From: %d ns \n', DispTimeAxisFrom);

DispTimeAxisTo = fread(fid, 1, 'int32');
fprintf(1,'        Time Axis To: %d ns \n', DispTimeAxisTo);

DispCountAxisFrom = fread(fid, 1, 'int32');
fprintf(1,'     Count Axis From: %d\n', DispCountAxisFrom); 

DispCountAxisTo = fread(fid, 1, 'int32');
fprintf(1,'       Count Axis To: %d\n', DispCountAxisTo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:8
DispCurveMapTo(i) = fread(fid, 1, 'int32');
DispCurveShow(i) = fread(fid, 1, 'int32');
fprintf(1,'---------------------------\n');
fprintf(1,'          Curve No %d\n', i-1);
fprintf(1,'               MapTo: %d\n', DispCurveMapTo(i));
fprintf(1,'                Show: %d\n', DispCurveShow(i));
fprintf(1,'---------------------------\n');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:3
ParamStart(i) = fread(fid, 1, 'float');
ParamStep(i) = fread(fid, 1, 'float');
ParamEnd(i) = fread(fid, 1, 'float');
fprintf(1,'---------------------------\n');
fprintf(1,'      Parameter No %d\n', i-1);
fprintf(1,'              Start: %d\n', ParamStart(i));
fprintf(1,'               Step: %d\n', ParamStep(i));
fprintf(1,'                End: %d\n', ParamEnd(i));
fprintf(1,'---------------------------\n');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RepeatMode = fread(fid, 1, 'int32');
fprintf(1,'        Repeat Mode: %d\n', RepeatMode);

RepeatsPerCurve = fread(fid, 1, 'int32');
fprintf(1,'     Repeat / Curve: %d\n', RepeatsPerCurve);

RepatTime = fread(fid, 1, 'int32');
fprintf(1,'        Repeat Time: %d\n', RepatTime);

RepeatWaitTime = fread(fid, 1, 'int32');
fprintf(1,'   Repeat Wait Time: %d\n', RepeatWaitTime);

ScriptName = setstr(fread(fid, 20, 'char'));
fprintf(1,'        Script Name: %s\n', ScriptName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          Header for each board
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1:NumberOfBoards
fprintf(1,'---------------------------\n');    
fprintf(1,'Board No %d\n', i-1);

HardwareIdent(:,i) = setstr(fread(fid, 16, 'char'));
fprintf(1,'Hardware Identifier: %s\n', HardwareIdent(:,i));

HardwareVersion(:,i) = setstr(fread(fid, 8, 'char'));
fprintf(1,'   Hardware Version: %s\n', HardwareVersion(:,i));    
    
HardwareSerial(i) = fread(fid, 1, 'int32');
fprintf(1,'   HW Serial Number: %d\n', HardwareSerial(i));

SyncDivider(i) = fread(fid, 1, 'int32');
fprintf(1,'       Sync divider: %d \n', SyncDivider(i));

CFDZeroCross0(i) = fread(fid, 1, 'int32');
fprintf(1,'    CFD0 zero cross: %d mV\n', CFDZeroCross0(i));

CFDLevel0(i) = fread(fid, 1, 'int32');
fprintf(1,'         CFD0 level: %d mV\n', CFDLevel0(i));

CFDZeroCross1(i) = fread(fid, 1, 'int32');
fprintf(1,'    CFD1 zero cross: %d mV\n', CFDZeroCross1(i));

CFDLevel1(i) = fread(fid, 1, 'int32');
fprintf(1,'         CFD1 level: %d mV\n', CFDLevel1(i));

Resolution(i) = fread(fid, 1, 'float');
fprintf(1,'         Resolution: %2.6g ns\n', Resolution(i));

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                Headers for each histogram (curve)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:NumberOfCurves

fprintf(1,'---------------------------\n');    
CurveIndex(i) = fread(fid, 1, 'int32');
fprintf(1,'      Curve Index: %d\n', CurveIndex(i));

TimeOfRecording(i) = fread(fid, 1, 'uint');

%  The PicoHarp software saves the time of recording
%  in a 32 bit serial time value as defined in all C libraries.
%  This equals the number of seconds elapsed since midnight
%  (00:00:00), January 1, 1970, coordinated universal time.
%  The conversion to normal date and time strings is everything but easy ....

TimeOfRecording(i) = TimeOfRecording(i)/24/60/60+25569+693960;
fprintf(1,'  Time of Recording: %s \n', datestr(TimeOfRecording(i),'dd-mmm-yyyy HH:MM:SS'));

HardwareIdent(:,i) = setstr(fread(fid, 16, 'char'));
fprintf(1,'Hardware Identifier: %s\n', HardwareIdent(:,i));
    
HardwareVersion(:,i) = setstr(fread(fid, 8, 'char'));
fprintf(1,'   Hardware Version: %s\n', HardwareVersion(:,i));    
    
HardwareSerial(i) = fread(fid, 1, 'int32');
fprintf(1,'   HW Serial Number: %d\n', HardwareSerial(i));

SyncDivider(i) = fread(fid, 1, 'int32');
fprintf(1,'       Sync divider: %d \n', SyncDivider(i));

CFDZeroCross0(i) = fread(fid, 1, 'int32');
fprintf(1,'    CFD0 zero cross: %d mV\n', CFDZeroCross0(i));

CFDLevel0(i) = fread(fid, 1, 'int32');
fprintf(1,'         CFD0 level: %d mV\n', CFDLevel0(i));

CFDZeroCross1(i) = fread(fid, 1, 'int32');
fprintf(1,'    CFD1 zero cross: %d mV\n', CFDZeroCross1(i));

CFDLevel1(i) = fread(fid, 1, 'int32');
fprintf(1,'         CFD1 level: %d mV\n', CFDLevel1(i));

Offset(i) = fread(fid, 1, 'int32');
fprintf(1,'             Offset: %d\n', Offset(i));

RoutingChannel(i) = fread(fid, 1, 'int32');
fprintf(1,'    Routing channel: %d\n', RoutingChannel(i));

ExtDevices(i) = fread(fid, 1, 'int32');
fprintf(1,'   External Devices: %d \n', ExtDevices(i));

MeasMode(i) = fread(fid, 1, 'int32');
fprintf(1,'   Measurement Mode: %d\n', MeasMode(i));

SubMode(i) = fread(fid, 1, 'int32');
fprintf(1,'           Sub-Mode: %d\n', SubMode(i));

P1(i) = fread(fid, 1, 'float');
fprintf(1,'                 P1: %d\n', P1(i));
P2(i) = fread(fid, 1, 'float');
fprintf(1,'                 P2: %d\n', P2(i));
P3(i) = fread(fid, 1, 'float');
fprintf(1,'                 P3: %d\n', P3(i));

RangeNo(i) = fread(fid, 1, 'int32');
fprintf(1,'          Range No.: %d\n', RangeNo(i));

Resolution(i) = fread(fid, 1, 'float');
fprintf(1,'         Resolution: %2.6g ns \n', Resolution(i));

Channels(i) = fread(fid, 1, 'int32');
fprintf(1,'           Channels: %d \n', Channels(i));

Tacq(i) = fread(fid, 1, 'int32');
fprintf(1,'   Acquisition Time: %d ms \n', Tacq(i));

StopAfter(i) = fread(fid, 1, 'int32');
fprintf(1,'         Stop After: %d ms \n', StopAfter(i));

StopReason(i) = fread(fid, 1, 'int32');
fprintf(1,'        Stop Reason: %d\n', StopReason(i));

InpRate0(i) = fread(fid, 1, 'int32');
fprintf(1,'       Input Rate 0: %d Hz\n', InpRate0(i));

InpRate1(i) = fread(fid, 1, 'int32');
fprintf(1,'       Input Rate 1: %d Hz\n', InpRate1(i));

HistCountRate(i) = fread(fid, 1, 'int32');
fprintf(1,'   Hist. Count Rate: %d cps\n', HistCountRate(i));

IntegralCount(i) = fread(fid, 1, 'int64');
fprintf(1,'     Integral Count: %d\n', IntegralCount(i));

Reserved(i) = fread(fid, 1, 'int32');
fprintf(1,'           Reserved: %d\n', Reserved(i));

DataOffset(i) = fread(fid, 1, 'int32');
fprintf(1,'         DataOffset: %d\n', DataOffset(i));

end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          Reads all histograms into one matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:NumberOfCurves
    fseek(fid,DataOffset(i),'bof');
    Counts(:,i) = fread(fid, Channels(i), 'uint32');
end;

Peak=max(Counts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          Summary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'\n');
fprintf(1,'\n');
fprintf(1,'=====================================================\n');
fprintf(1,'                     SUMMARY                         \n');
fprintf(1,'=====================================================\n');
fprintf(1,' Curve    Channel     Number of    Peak     Integral \n');
fprintf(1,' index   resolution   channels     count     count   \n');
fprintf(1,'=====================================================\n');

for i = 1:NumberOfCurves
fprintf(1,'  %3i       %2.6g  %10i  %10i  %10i\n', CurveIndex(i),Resolution(i), Channels(i), Peak(i), IntegralCount(i));   
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          This is a simple display of the histogram(s)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(1);
semilogy(Counts);
% axis([0 max(max(Channels)) 1 10*max(max(Counts))]);
xlabel('Channel #');
ylabel('Counts');

if NumberOfCurves<21
   legend(num2str((1:NumberOfCurves)'),0);
end;

fclose(fid);