% Demo accessing TimeHarp 200 continuous mode data files (*.thc) from MATLAB
% TimeHarp 200, Software version 2.0 through 5.1
% Tested with Matlab 4 and 5.
% Peter Kapusta, PicoQuant GmbH 2003
% This is demo code. Use at your own risk. No warranties.

clear all;
clc;
fprintf(1,'\n');

[filename, pathname]=uigetfile('*.thc', 'Continuous mode data:', 0, 0);
fid=fopen(filename);

%
% Read the TxtHdr Record
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
% Read the BinHdr Record
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
% Read the BoardHdr Record
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

end;



%
% Channels
%


FirstChannel = fread(fid, 1, 'int32');
fprintf(1,'First Channel: %d\n', FirstChannel);

LastChannel = fread(fid, 1, 'int32');
fprintf(1,'Last Channel: %d\n', LastChannel);

%
% Curves
%

for i = 1:Bin_NumberOfCurves
	for j = 1:Bin_NumberOfChannels
	Counts(i,j) = fread(fid, 1, 'uint16');
	end;
end;


figure(1);

%semilogy(1:Bin_NumberOfChannels,Counts);
semilogy(FirstChannel:LastChannel,Counts);
xlabel('Channel #');
ylabel('Counts');


figure(2);
[X,Y]=meshgrid(FirstChannel:LastChannel,1:Bin_NumberOfCurves);
shading flat;
mesh(X,Y,log10(Counts));
ylabel('Histogram #'),xlabel('Channel #'),zlabel('Log Counts');

fclose(fid);