% Demo accessing TimeHarp interactive mode data files (*.thd) from MATLAB
% TimeHarp 200, Software version 2.0 through 5.2
% Tested with Matlab 5 through 6.
% Peter Kapusta, PicoQuant GmbH 2003
% This is demo code. Use at your own risk. No warranties.
% Make sure you have enough memory when loading large files!

clear all;
clc;
fprintf(1,'\n');

[filename, pathname]=uigetfile('*.thd', 'Interactive mode data:', 0, 0)

fid=fopen([pathname filename]);
%
% Read the TxtHdr
%

Txt_Ident = setstr(fread(fid, 16, 'char'));
fprintf(1,'Ident: %s\n', Txt_Ident);

Txt_SoftwareVersion = setstr(fread(fid, 6, 'char'));
fprintf(1,'Software version: %s\n', Txt_SoftwareVersion);

Txt_HardWareVersion = setstr(fread(fid, 6, 'char'));
fprintf(1,'Hardware version: %s\n', Txt_HardWareVersion);

Txt_FileTime = setstr(fread(fid, 18, 'char'));
fprintf(1,'File creation: %s\n', Txt_FileTime);

Txt_CRLF = setstr(fread(fid, 2, 'char'));

Txt_Comment = setstr(fread(fid, 256, 'char'));
fprintf(1,'Comment: %s\n', Txt_Comment);


%
% Read the BinHdr
%


Bin_NumberOfChannels = fread(fid, 1, 'int32');
fprintf(1,'Number of Channels: %d\n', Bin_NumberOfChannels);

Bin_NumberOfCurves = fread(fid, 1, 'int32');
fprintf(1,'Number of Curves: %d\n', Bin_NumberOfCurves);

Bin_BitsPerChannel = fread(fid, 1, 'int32');
fprintf(1,'Bits / Channel: %d\n', Bin_BitsPerChannel);

Bin_RoutingChannels = fread(fid, 1, 'int32');
fprintf(1,'Routing Channels: %d\n', Bin_RoutingChannels);

Bin_NumberOfBoards = fread(fid, 1, 'int32');
fprintf(1,'Number of Boards: %d\n', Bin_NumberOfBoards);

Bin_ActiveCurve = fread(fid, 1, 'int32');
fprintf(1,'Active Curve: %d\n', Bin_ActiveCurve);

Bin_MeasurementMode = fread(fid, 1, 'int32');
fprintf(1,'Measurement Mode: %d\n', Bin_MeasurementMode);

Bin_SubMode = fread(fid, 1, 'int32');
fprintf(1,'SubMode: %d\n', Bin_SubMode);

Bin_RangeNo = fread(fid, 1, 'int32');
fprintf(1,'Range No.: %d\n', Bin_RangeNo);

Bin_Offset = fread(fid, 1, 'int32');
fprintf(1,'Offset: %d ns \n', Bin_Offset);

Bin_AcquisitionTime = fread(fid, 1, 'int32');
fprintf(1,'Acquisition Time: %d ms \n', Bin_AcquisitionTime);

Bin_StopAt = fread(fid, 1, 'int32');
fprintf(1,'Stop at: %d counts \n', Bin_StopAt);

Bin_StopOnOvfl = fread(fid, 1, 'int32');
fprintf(1,'Stop on Overflow: %d\n', Bin_StopOnOvfl);

Bin_Restart = fread(fid, 1, 'int32');
fprintf(1,'Restart: %d\n', Bin_Restart);

Bin_DispLinLog = fread(fid, 1, 'int32');
fprintf(1,'Display Lin/Log: %d\n', Bin_DispLinLog);

Bin_DispTimeAxisFrom = fread(fid, 1, 'int32');
fprintf(1,'Display Time Axis From: %d ns \n', Bin_DispTimeAxisFrom);

Bin_DispTimeAxisTo = fread(fid, 1, 'int32');
fprintf(1,'Display Time Axis To: %d ns \n', Bin_DispTimeAxisTo);

Bin_DispCountAxisFrom = fread(fid, 1, 'int32');
fprintf(1,'Display Count Axis From: %d\n', Bin_DispCountAxisFrom); 

Bin_DispCountAxisTo = fread(fid, 1, 'int32');
fprintf(1,'Display Count Axis To: %d\n', Bin_DispCountAxisTo);

for i = 1:8
Bin_DispCurveMapTo(i) = fread(fid, 1, 'int32');
Bin_DispCurveShow(i) = fread(fid, 1, 'int32');
end;

Bin_DispCurve=0:7
Bin_DispCurveMapTo
Bin_DispCurveShow

for i = 1:3
Bin_ParamStart(i) = fread(fid, 1, 'float');
Bin_ParamStep(i) = fread(fid, 1, 'float');
Bin_ParamEnd(i) = fread(fid, 1, 'float');
end;

Bin_ParamStart
Bin_ParamStep
Bin_ParamEnd

Bin_RepeatMode = fread(fid, 1, 'int32');
fprintf(1,'Repeat Mode: %d\n', Bin_RepeatMode);

Bin_RepeatsPerCurve = fread(fid, 1, 'int32');
fprintf(1,'Repeat / Curve: %d\n', Bin_RepeatsPerCurve);

Bin_RepatTime = fread(fid, 1, 'int32');
fprintf(1,'Repeat Time: %d\n', Bin_RepatTime);

Bin_RepeatWait = fread(fid, 1, 'int32');
fprintf(1,'Repeat Wait Time: %d\n', Bin_RepeatWait);

Bin_ScriptName = setstr(fread(fid, 20, 'char'));
fprintf(1,'Script Name: %s\n', Bin_ScriptName);


%
% Read the BoardHdr
%


Board_BoardSerial = fread(fid, 1, 'int32');
fprintf(1,'Board Serial Number: %d\n', Board_BoardSerial);

Board_CFDZeroCross = fread(fid, 1, 'int32');
fprintf(1,'CFD Zero Cross: %d mV\n', Board_CFDZeroCross);

Board_CFDDiscriminatorMin = fread(fid, 1, 'int32');
fprintf(1,'CFD Discriminator Min.: %d mV\n', Board_CFDDiscriminatorMin);

Board_SYNCLevel = fread(fid, 1, 'int32');
fprintf(1,'SYNC Level: %d mV\n', Board_SYNCLevel);

Board_CurveOffset = fread(fid, 1, 'int32');
fprintf(1,'Curve Offset: %d\n', Board_CurveOffset);

Board_Resolution = fread(fid, 1, 'float');
fprintf(1,'Resolution: %5.6f ns\n', Board_Resolution);


%
% Read the Stored Curve Records
%


for i = 1:Bin_NumberOfCurves

fprintf(1,'\n');

CurveIndex(i) = fread(fid, 1, 'int32');
fprintf(1,'Curve Index: %d\n', CurveIndex(i));

TimeOfRecording(i) = fread(fid, 1, 'int');

%  The TimeHarp software saves the time of recording
%  in a 32 bit time value as defined in all C libraries.
%  This equals the number of seconds elapsed since midnight
%  (00:00:00), January 1, 1970, coordinated universal time.

%  The following (optional) code line converts the 32 bit time value
%  to a value of TDateTime type, used for example in Delphi.
 
TimeOfRecording(i) = TimeOfRecording(i)/24/60/60+25569;

%  The integral part of the new value is the number of days that
%  have passed since December 30, 1899. The fractional part of 
%  the number is a fraction of a 24 hour day that has elapsed.

fprintf(1,'Time of Recording: %d \n', TimeOfRecording(i));

BoardSerial(i) = fread(fid, 1, 'int32');
fprintf(1,'Board Serial Number: %d\n', BoardSerial(i));

CFDZeroCross(i) = fread(fid, 1, 'int32');
fprintf(1,'CFD Zero Cross: %d mV\n', CFDZeroCross(i));

CFDDiscriminatorMin(i) = fread(fid, 1, 'int32');
fprintf(1,'CFD Discriminator Min.: %d mV\n', CFDDiscriminatorMin(i));

SYNCLevel(i) = fread(fid, 1, 'int32');
fprintf(1,'SYNC Level: %d mV\n', SYNCLevel(i));

CurveOffset(i) = fread(fid, 1, 'int32');
fprintf(1,'Curve Offset: %d\n', CurveOffset(i));

RoutingChannel(i) = fread(fid, 1, 'int32');
fprintf(1,'Routing Channel: %d\n', RoutingChannel(i));

SubMode(i) = fread(fid, 1, 'int32');
fprintf(1,'SubMode: %d\n', SubMode(i));

MeasMode(i) = fread(fid, 1, 'int32');
fprintf(1,'Measurement Mode: %d\n', MeasMode(i));

P1(i) = fread(fid, 1, 'float');
fprintf(1,'P1: %d\n', P1(i));
P2(i) = fread(fid, 1, 'float');
fprintf(1,'P2: %d\n', P2(i));
P3(i) = fread(fid, 1, 'float');
fprintf(1,'P3: %d\n', P3(i));

RangeNo(i) = fread(fid, 1, 'int32');
fprintf(1,'Range No.: %d\n', RangeNo(i));

Offset(i) = fread(fid, 1, 'int32');
fprintf(1,'Offset: %d ns \n', Offset(i));

AcquisitionTime(i) = fread(fid, 1, 'int32');
fprintf(1,'Acquisition Time: %d ms \n', AcquisitionTime(i));

StopAfter(i) = fread(fid, 1, 'int32');
fprintf(1,'Stop After: %d ms \n', StopAfter(i));

StopReason(i) = fread(fid, 1, 'int32');
fprintf(1,'Stop Reason: %d\n', StopReason(i));

SyncRate(i) = fread(fid, 1, 'int32');
fprintf(1,'SYNC Rate: %d Hz\n', SyncRate(i));

CFDCountRate(i) = fread(fid, 1, 'int32');
fprintf(1,'CFD Count Rate: %d cps\n', CFDCountRate(i));

TDCCountRate(i) = fread(fid, 1, 'int32');
fprintf(1,'TDC Count Rate: %d cps\n', TDCCountRate(i));

IntegralCount(i) = fread(fid, 1, 'int32');
fprintf(1,'Integral Count: %d\n', IntegralCount(i));

Resolution(i) = fread(fid, 1, 'float');
fprintf(1,'Resolution: %5.6f ns\n', Resolution(i));

Extdev(i) = fread(fid, 1, 'int32');
fprintf(1,'External Device: %d\n', Extdev(i));

Reserved(i) = fread(fid, 1, 'int32');
fprintf(1,'Reserved: %d\n', Reserved(i));

for j = 1:Bin_NumberOfChannels
   Counts(i,j) = fread(fid, 1, 'uint32');
end;

end;


figure(1);
semilogy(1:Bin_NumberOfChannels,Counts);
xlabel('Channel #');
ylabel('Counts');



fclose(fid);