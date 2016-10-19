function [x, head] = p3rRead(name,cnt);

% [x, head] = p3rRead(name,cnts);

if ~(name(end-2:end)=='p3r')
    disp('Not a p3r-file!');
    return
else
    fin = fopen(name,'r');
    if (fin==-1)
        errordlg('Cannot open specified file. Please try again.');
    else

        head = struct('Ident', setstr(fread(fin, 16, 'char')));
        head = setfield(head, 'FormatVersion', setstr(fread(fin, 6, 'char')));
        head = setfield(head, 'CreatorName', setstr(fread(fin, 18, 'char')));
        head = setfield(head, 'CreatorVersion', setstr(fread(fin, 12, 'char')));
        head = setfield(head, 'FileTime', setstr(fread(fin, 18, 'char')));
        head = setfield(head, 'CRLF', setstr(fread(fin, 2, 'char')));
        head = setfield(head, 'CommentField', setstr(fread(fin, 256, 'char')));

        head = setfield(head, 'Curves', fread(fin, 1, 'int32'));
        head = setfield(head, 'BitsPerRecord', fread(fin, 1, 'int32'));
        head = setfield(head, 'RoutingChannels', fread(fin, 1, 'int32'));
        head = setfield(head, 'NumberOfBoards', fread(fin, 1, 'int32'));
        head = setfield(head, 'ActiveCurve', fread(fin, 1, 'int32'));
        head = setfield(head, 'MeasurementMode', fread(fin, 1, 'int32'));
        head = setfield(head, 'SubMode', fread(fin, 1, 'int32'));
        head = setfield(head, 'RangeNo', fread(fin, 1, 'int32'));
        head = setfield(head, 'Offset', fread(fin, 1, 'int32'));
        head = setfield(head, 'AcquisitionTime', fread(fin, 1, 'int32'));
        head = setfield(head, 'StopAt', fread(fin, 1, 'int32'));
        head = setfield(head, 'StopOnOvfl', fread(fin, 1, 'int32'));
        head = setfield(head, 'Restart', fread(fin, 1, 'int32'));
        head = setfield(head, 'DispLinLog', fread(fin, 1, 'int32'));
        head = setfield(head, 'DispTimeFrom', fread(fin, 1, 'int32'));
        head = setfield(head, 'DispTimeTo', fread(fin, 1, 'int32'));
        head = setfield(head, 'DispCountFrom', fread(fin, 1, 'int32'));
        head = setfield(head, 'DispCountTo', fread(fin, 1, 'int32'));

        for i = 1:8
            DispCurveMapTo(i) = fread(fin, 1, 'int32');
            DispCurveShow(i) = fread(fin, 1, 'int32');
        end;
        DispCurve=0:7;

        head = setfield(head, 'DispCurve', DispCurve);
        head = setfield(head, 'DispCurveMapTo', DispCurveMapTo);
        head = setfield(head, 'DispCurveShow', DispCurveShow);

        for i = 1:3
            ParamStart(i) = fread(fin, 1, 'float');
            ParamStep(i) = fread(fin, 1, 'float');
            ParamEnd(i) = fread(fin, 1, 'float');
        end;

        head = setfield(head, 'ParamStart', ParamStart);
        head = setfield(head, 'ParamStep', ParamStep);
        head = setfield(head, 'ParamEnd', ParamEnd);

        head = setfield(head, 'RepeatMode', fread(fin, 1, 'int32'));
        head = setfield(head, 'RepeatsPerCurve', fread(fin, 1, 'int32'));
        head = setfield(head, 'RepeatTime', fread(fin, 1, 'int32'));
        head = setfield(head, 'RepeatWait', fread(fin, 1, 'int32'));
        head = setfield(head, 'ScriptName', setstr(fread(fin, 20, 'char')));

        head = setfield(head, 'HardwareIdent', setstr(fread(fin, 16, 'char')));
        head = setfield(head, 'HardwareVersion', setstr(fread(fin, 8, 'char')));
        head = setfield(head, 'HardwareSerial', fread(fin, 1, 'int32'));
        head = setfield(head, 'SyncDivider', fread(fin, 1, 'int32'));
        head = setfield(head, 'CFDZeroCross0', fread(fin, 1, 'int32'));
        head = setfield(head, 'CFDLevel0', fread(fin, 1, 'int32'));
        head = setfield(head, 'CFDZeroCross1', fread(fin, 1, 'int32'));
        head = setfield(head, 'CFDLevel1', fread(fin, 1, 'int32'));
        head = setfield(head, 'Resolution', fread(fin, 1, 'float'));

        head = setfield(head, 'ExtDevices', fread(fin, 1, 'int32'));
        head = setfield(head, 'Reserved1', fread(fin, 1, 'int32'));
        head = setfield(head, 'Reserved2', fread(fin, 1, 'int32'));
        head = setfield(head, 'CntRate0', fread(fin, 1, 'int32'));
        head = setfield(head, 'CntRate1', fread(fin, 1, 'int32'));
        head = setfield(head, 'StopAfter', fread(fin, 1, 'int32'));
        head = setfield(head, 'StopReason', fread(fin, 1, 'int32'));
        head = setfield(head, 'Records', fread(fin, 1, 'uint32'));
        head = setfield(head, 'ImgHdrSize', fread(fin, 1, 'int32'));

        %Special header for imaging - currently not used
        %SpecHeader = fread(fin, ImgHdrSize, 'int32');


        %
        %  This reads the T3 mode event records
        %

        ofltime = 0;
        dlen = 0;
        WRAPAROUND=65536;

        syncperiod = 1E9/head.CntRate0;   % in nanoseconds

        if nargin<2
            cnt = [1 inf];
        end

        fseek(fin, 4*(cnt(1)-1), 0);
        x = fread(fin, cnt(2), 'ubit32');

        chan = floor(x/2^28);
        dtime = floor(x/2^16 - chan*2^12);
        x = x - 2^16*floor(x/2^16);
        x = x + WRAPAROUND*cumsum(chan==15 & dtime ==0);
        x(~(chan==1)) = [];
        dtime(~(chan==1)) = [];
        chan(~(chan==1)) = [];
        x = x*syncperiod + dtime*head.Resolution;

        fclose(fin);
    end
end
