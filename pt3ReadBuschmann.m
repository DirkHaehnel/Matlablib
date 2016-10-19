function [head, y, tcspc, z, num, overcount] = pt3ReadBuschmann(name, start_cnt, start_time, max_time, sync);

% [head, y, tcspc, num, overcount] = pt3Read(name,[start_cnt, start_time, max_time, sync);

if ~(name(end-2:end)=='pt3')
    disp('Not a pt3-file!');
    return
else
    fin = fopen(name,'r');
    if (fin==-1)
        errordlg('Cannot open specified file. Please try again.');
    else
        head = struct('Ident', char(fread(fin, 16, 'char')'));
        head = setfield(head, 'FormatVersion', char(fread(fin, 6, 'char')'));
        head = setfield(head, 'CreatorName', char(fread(fin, 18, 'char')'));
        head = setfield(head, 'CreatorVersion', char(fread(fin, 12, 'char')'));
        head = setfield(head, 'FileTime', char(fread(fin, 18, 'char')));
        fread(fin, 2, 'char');
        head = setfield(head, 'Comment', char(fread(fin, 256, 'char')'));

        head = setfield(head, 'NCurves', fread(fin, 1, 'int32'));
        head = setfield(head, 'NChannels', 4096);
        head = setfield(head, 'BitsPerRecord', fread(fin, 1, 'int32'));
        head = setfield(head, 'RoutingChannels', fread(fin, 1, 'int32'));
        head = setfield(head, 'NumberOfBoards', fread(fin, 1, 'int32'));
        head = setfield(head, 'ActiveCurve', fread(fin, 1, 'int32'));
        head = setfield(head, 'MeasMode', fread(fin, 1, 'int32'));
        head = setfield(head, 'SubMode', fread(fin, 1, 'int32'));
        head = setfield(head, 'RangeNo', fread(fin, 1, 'int32'));
        head = setfield(head, 'Offset', fread(fin, 1, 'int32'));
        head = setfield(head, 'TAcq', fread(fin, 1, 'int32'));
        head = setfield(head, 'StopAt', fread(fin, 1, 'int32'));
        head = setfield(head, 'StopOnOvfl', fread(fin, 1, 'int32'));
        head = setfield(head, 'Restart', fread(fin, 1, 'int32'));
        head = setfield(head, 'LinLog', fread(fin, 1, 'int32'));
        head = setfield(head, 'MinAx', fread(fin, 1, 'int32'));
        head = setfield(head, 'MaxAx', fread(fin, 1, 'int32'));
        head = setfield(head, 'MinAxCnt', fread(fin, 1, 'int32'));
        head = setfield(head, 'MaxAxCnt', fread(fin, 1, 'int32'));
        head = setfield(head, 'DispCurves', fread(fin, 16, 'int32'));
        head = setfield(head, 'Params', fread(fin, 9, 'int32'));
        head = setfield(head, 'RepeatMode', fread(fin, 1, 'int32'));
        head = setfield(head, 'RepeatsPerCurve', fread(fin, 1, 'int32'));
        head = setfield(head, 'RepeatTime', fread(fin, 1, 'int32'));
        head = setfield(head, 'RepeatWaitTime', fread(fin, 1, 'int32'));
        head = setfield(head, 'ScriptName', fread(fin, 20, 'char'));
        head = setfield(head, 'BordIdent',  fread(fin, 16, 'char'));
        head = setfield(head, 'BordVersion', fread(fin, 8, 'char'));
        head = setfield(head, 'BordSerial', fread(fin, 1, 'int32'));
        head = setfield(head, 'SyncDiv', fread(fin, 1, 'int32'));
        head = setfield(head, 'CFDZeroCross0', fread(fin, 1, 'int32'));
        head = setfield(head, 'CFDLevel0', fread(fin, 1, 'int32'));
        head = setfield(head, 'CFDZeroCross1', fread(fin, 1, 'int32'));
        head = setfield(head, 'CFDLevel1', fread(fin, 1, 'int32'));
        head = setfield(head, 'Resolution', fread(fin, 1, 'float'));

        head = setfield(head, 'ExternalDev', fread(fin, 1, 'int32'));
        head = setfield(head, 'Reserved1', fread(fin, 1, 'int32'));
        head = setfield(head, 'Reserved2', fread(fin, 1, 'int32'));
        head = setfield(head, 'CntRate0', fread(fin, 1, 'int32'));
        head = setfield(head, 'CntRate1', fread(fin, 1, 'int32'));
        head = setfield(head, 'StopAfter', fread(fin, 1, 'int32'));
        head = setfield(head, 'StopReason', fread(fin, 1, 'int32'));
        head = setfield(head, 'NCounts', fread(fin, 1, 'uint32'));

        SpecHeaderLength = fread(fin, 1, 'int32');

        if SpecHeaderLength>0
            tmp = fread(fin, 2, 'uint32');
            head = setfield(head, 'ScanIdent', tmp(2));
            if tmp(2)==1 % PI E710 Scan Controller
                head = setfield(head, 'ScanTimePerPix', fread(fin, 1, 'int32'));
                fread(fin, 1, 'int32');
                head = setfield(head, 'ScanPattern', fread(fin, 1, 'int32'));
                fread(fin, 1, 'int32');
                head = setfield(head, 'ScanStartX', fread(fin, 1, 'float'));
                head = setfield(head, 'ScanStartY', fread(fin, 1, 'float'));
                head = setfield(head, 'ScanWidthX', fread(fin, 1, 'int32'));
                head = setfield(head, 'ScanWidthY', fread(fin, 1, 'int32'));
                head = setfield(head, 'ScanResolution', fread(fin, 1, 'float'));
                head = setfield(head, 'ScanTStartTo', fread(fin, 1, 'float'));
                head = setfield(head, 'ScanTStopTo', fread(fin, 1, 'float'));
                head = setfield(head, 'ScanTStartFro', fread(fin, 1, 'float'));
                head = setfield(head, 'ScanTStopFro', fread(fin, 1, 'float'));
            end
            if tmp(2)==2 % SCX 200 Scan Controller
                head = setfield(head, 'ScanTimePerPix', fread(fin, 1, 'int32'));
                head = setfield(head, 'ScanPause', fread(fin, 1, 'int32'));
                head = setfield(head, 'ScanPattern', fread(fin, 1, 'int32'));
                fread(fin, 1, 'int32');
                head = setfield(head, 'ScanStartX', fread(fin, 1, 'float'));
                head = setfield(head, 'ScanStartY', fread(fin, 1, 'float'));
                head = setfield(head, 'ScanWidthX', fread(fin, 1, 'int32'));
                head = setfield(head, 'ScanWidthY', fread(fin, 1, 'int32'));
                fread(fin, 1, 'int32');
            end
        end

        if nargin>1
            head.CntRate0 = sync;
            
            fseek(fin, 4*(start_cnt-1), 0);
            
            ind       = [];
            valid     = [];
            overflow  = [];
            tcspc     = [];
            y         = [];
            z         = [];
            num       = 0;
            overcount = start_time;
            sync = 1e9/head.CntRate0;
            
            while (sum(ind)==0)&(~feof(fin))

                [tmpy num] = fread(fin, 5E5, 'uint32');
                
                chan    = floor(tmpy/2^28);
                tmpt    = floor(tmpy/2^16 - chan*2^12);
                numsync = tmpy - 2^16*floor(tmpy/2^16);   
                
                valid     = [valid; (chan==0)];
                overflow  = [overflow; (chan==15&tmpt==0)];
                
                cs_ovl              = start_time + 2^16*cumsum(overflow);
                cs_ovl(1:length(y)) = [];
                
                y     = [y; sync*(numsync + cs_ovl)];
                tcspc = [tcspc; tmpt*head.Resolution];            % convert to ns;
                z     = [z; chan];
                               
                ind    = (y>=max_time);
            end;
            
            if (num>0)
                y(ind)        = [];
                tcspc(ind)    = [];
                z(ind)        = [];
                valid(ind)    = [];
                overflow(ind) = [];

                num       = length(y);
                cs_ovl    = start_time + 2^16*cumsum(overflow);
                if length(cs_ovl)>0
                    overcount = cs_ovl(end);
                end;

                y(valid==1)     = [];
                tcspc(valid==1) = [];
                z(valid==1)     = [];
            end;
            
        end
        fclose(fin);
    end
end
