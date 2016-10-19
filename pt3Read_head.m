function [head, HeaderLength] = pt3Read_head(name)

% [head, HeaderLength] = pt3Read_head(name);

head = [];

if strcmp(name(end-3:end),'.pt3')
    fin = fopen(name,'r');
    if (fin==-1)
        errordlg('Cannot open specified file. Please try again.');
    else
        HeaderLength = fread(fin, 1, 'int32');
        
        if HeaderLength > 0
            fseek(fin,0,'bof');
        end;
        
        head.Ident = fread(fin, 16, '*char')';
        head.FormatVersion = fread(fin, 6, '*char')';
        head.CreatorName = fread(fin, 18, '*char')';
        head.CreatorVersion = fread(fin, 12, '*char')';
        head.FileTime = fread(fin, 18, '*char')';
        fread(fin, 2, 'char');
        head.Comment = fread(fin, 256, '*char')';

        head.NCurves = fread(fin, 1, 'int32');
        head.NChannels = 4096;
        head.BitsPerRecord = fread(fin, 1, 'int32');
        head.RoutingChannels = fread(fin, 1, 'int32');
        head.NumberOfBoards = fread(fin, 1, 'int32');
        head.ActiveCurve = fread(fin, 1, 'int32');
        head.MeasMode = fread(fin, 1, 'int32');
        head.SubMode = fread(fin, 1, 'int32');
        head.RangeNo = fread(fin, 1, 'int32');
        head.Offset = fread(fin, 1, 'int32');
        head.TAcq = fread(fin, 1, 'int32');
        head.StopAt = fread(fin, 1, 'int32');
        head.StopOnOvfl = fread(fin, 1, 'int32');
        head.Restart = fread(fin, 1, 'int32');
        head.LinLog = fread(fin, 1, 'int32');
        head.MinAx = fread(fin, 1, 'int32');
        head.MaxAx = fread(fin, 1, 'int32');
        head.MinAxCnt = fread(fin, 1, 'int32');
        head.MaxAxCnt = fread(fin, 1, 'int32');
        head.DispCurves = fread(fin, 16, 'int32');
        head.Params = fread(fin, 9, 'int32');
        head.RepeatMode = fread(fin, 1, 'int32');
        head.RepeatsPerCurve = fread(fin, 1, 'int32');
        head.RepeatTime = fread(fin, 1, 'int32');
        head.RepeatWaitTime = fread(fin, 1, 'int32');
        
        head.ScriptName = fread(fin, 20, '*char')';
        
        for n=(1:head.NumberOfBoards)
            head.BoardIdent(:,n)   = fread(fin, 16, '*char')';
            head.BoardVersion(:,n) = fread(fin, 8, '*char')';
            head.BoardSerial(n)    = fread(fin, 1, 'int32');
            head.SyncDiv(n)        = fread(fin, 1, 'int32');
            head.CFDZeroCross0(n)  = fread(fin, 1, 'int32');
            head.CFDLevel0(n)      = fread(fin, 1, 'int32');
            head.CFDZeroCross1(n)  = fread(fin, 1, 'int32');
            head.CFDLevel1(n)      = fread(fin, 1, 'int32');
            head.Resolution(n)     = fread(fin, 1, 'float');
            if strcmp(head.FormatVersion(1:3),'2.0')
                head.RouterModel(n)        = fread(fin, 1, 'int32');
                head.RouterEnabled(n)      = fread(fin, 1, 'int32');
                head.RtCh1_InputType(n)    = fread(fin, 1, 'int32');
                head.RtCh1_InputLevel(n)   = fread(fin, 1, 'int32');
                head.RtCh1_InputEdge(n)    = fread(fin, 1, 'int32');
                head.RtCh1_CFDPresent(n)   = fread(fin, 1, 'int32');
                head.RtCh1_CFDLevel(n)     = fread(fin, 1, 'int32');
                head.RtCh1_CFDZeroCross(n) = fread(fin, 1, 'int32');
                head.RtCh2_InputType(n)    = fread(fin, 1, 'int32');
                head.RtCh2_InputLevel(n)   = fread(fin, 1, 'int32');
                head.RtCh2_InputEdge(n)    = fread(fin, 1, 'int32');
                head.RtCh2_CFDPresent(n)   = fread(fin, 1, 'int32');
                head.RtCh2_CFDLevel(n)     = fread(fin, 1, 'int32');
                head.RtCh2_CFDZeroCross(n) = fread(fin, 1, 'int32');
                head.RtCh3_InputType(n)    = fread(fin, 1, 'int32');
                head.RtCh3_InputLevel(n)   = fread(fin, 1, 'int32');
                head.RtCh3_InputEdge(n)    = fread(fin, 1, 'int32');
                head.RtCh3_CFDPresent(n)   = fread(fin, 1, 'int32');
                head.RtCh3_CFDLevel(n)     = fread(fin, 1, 'int32');
                head.RtCh3_CFDZeroCross(n) = fread(fin, 1, 'int32');
                head.RtCh4_InputType(n)    = fread(fin, 1, 'int32');
                head.RtCh4_InputLevel(n)   = fread(fin, 1, 'int32');
                head.RtCh4_InputEdge(n)    = fread(fin, 1, 'int32');
                head.RtCh4_CFDPresent(n)   = fread(fin, 1, 'int32');
                head.RtCh4_CFDLevel(n)     = fread(fin, 1, 'int32');
                head.RtCh4_CFDZeroCross(n) = fread(fin, 1, 'int32');
            end
        end;
        head.ExternalDev = fread(fin, 1, 'int32');
        head.Reserved1 = fread(fin, 1, 'int32');
        head.Reserved2 = fread(fin, 1, 'int32');
        head.CntRate0 = fread(fin, 1, 'int32');
        head.CntRate1 = fread(fin, 1, 'int32');
        head.StopAfter = fread(fin, 1, 'int32');
        head.StopReason = fread(fin, 1, 'int32');
        head.NCounts = fread(fin, 1, 'uint32');

        head.SpecHeaderLength = fread(fin, 1, 'int32');

        if head.SpecHeaderLength>0
            head.ImgHdr.Dimensions = fread(fin, 1, 'int32');
            head.ImgHdr.Ident = fread(fin, 1, 'int32');
            if head.ImgHdr.Ident==1 % PI E710 Scan Controller
                head.ImgHdr.ScanTimePerPix = fread(fin, 1, 'int32');
                fread(fin, 1, 'int32');
                head.ImgHdr.ScanPattern = fread(fin, 1, 'int32');
                fread(fin, 1, 'int32');
                head.ImgHdr.ScanStartX = fread(fin, 1, 'float');
                head.ImgHdr.ScanStartY = fread(fin, 1, 'float');
                head.ImgHdr.ScanWidthX = fread(fin, 1, 'int32');
                head.ImgHdr.ScanWidthY = fread(fin, 1, 'int32');
                head.ImgHdr.ScanResolution = fread(fin, 1, 'float');
                head.ImgHdr.ScanTStartTo = fread(fin, 1, 'float');
                head.ImgHdr.ScanTStopTo = fread(fin, 1, 'float');
                head.ImgHdr.ScanTStartFro = fread(fin, 1, 'float');
                head.ImgHdr.ScanTStopFro = fread(fin, 1, 'float');
                if head.SpecHeaderLength==16
                    head.ImgHdr.ScanOffset = fread(fin, 1, 'uint32');
                else
                    head.SpecHeaderLength = 16;
                    head.ImgHdr.ScanOffset = 0;
                end;
            end
            if head.ImgHdr.Ident==2 % SCX 200 Scan Controller
                head.ImgHdr.ScanTimePerPix = fread(fin, 1, 'int32');
                head.ImgHdr.ScanPause = fread(fin, 1, 'int32');
                head.ImgHdr.ScanPattern = fread(fin, 1, 'int32');
                fread(fin, 1, 'int32');
                head.ImgHdr.ScanStartX = fread(fin, 1, 'float');
                head.ImgHdr.ScanStartY = fread(fin, 1, 'float');
                head.ImgHdr.ScanWidthX = fread(fin, 1, 'int32');
                head.ImgHdr.ScanWidthY = fread(fin, 1, 'int32');
                if head.SpecHeaderLength==11
                    head.ImgHdr.ScanOffset = fread(fin, 1, 'uint32');
                else
                    head.SpecHeaderLength = 11;
                    head.ImgHdr.ScanOffset = 0;
                    fread(fin, 1, 'int32');
                end;
            end
            if head.ImgHdr.Ident==3 % LSM
                head.ImgHdr.Frame = fread(fin, 1, 'int32');
                head.ImgHdr.LineStart  = fread(fin, 1, 'int32');
                head.ImgHdr.LineStop   = fread(fin, 1, 'int32');
                head.ImgHdr.Pattern    = fread(fin, 1, 'int32');
                head.ImgHdr.ScanWidthX = fread(fin, 1, 'int32');
                head.ImgHdr.ScanWidthY = fread(fin, 1, 'int32');
                for n = 9:head.SpecHeaderLength
                    tmp = fread(fin, 1, 'uint32');
                end
            end
        end
        HeaderLength = ftell(fin);
    end;
    fclose(fin);
else
    disp('Not a pt3-file!');
end;



