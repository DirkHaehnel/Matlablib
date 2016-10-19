function [y, tcspc, flag, num, overcount, head] = pt3Read_s(name, start_cnt, start_time, sync)

%        [y, tcspc, flag, num, overcount, head] = pt3Read_s(name, start_cnt, start_time, sync)


if strcmp(name(end-3:end),'.pt3')
    hlength = 0;

    if nargin>1

        if (nargin<4)||(nargout==6) 
            [head, hlength] = pt3Read_head(name);
            sync            = 1e9/head.CntRate0;
        end;

        if nargin<3
            start_time = 0;
        end;
        
        if start_cnt<1
            start_cnt = 1;
        end;

    end;

    fin = fopen(name,'r');
    if (fin==-1)
        errordlg('Cannot open specified file. Please try again.');
    else
        if (hlength==0)
            hlength = fread(fin, 1, 'int32');
            if hlength > 0
                fclose(fin);
                pt3_headerupdate(name);
                fin = fopen(name,'r');
                hlength = fread(fin, 1, 'int32');
            end;
            hlength = abs(hlength);
        end

        fseek(fin, hlength+4*(start_cnt-1), 'bof');

        valid     = [];
        overflow  = [];
        y         = [];
        tcspc     = [];
        flag      = [];

        num       = 0;
        overcount = start_time;

        mark = [];

        while (length(mark)<5)

            [tmpy num] = fread(fin, 5E5, 'uint32');

            tmpy = uint32(tmpy);

            chan    = bitshift(tmpy, -28);
            tmpt    = double(bitshift(bitand(tmpy,268369920),-16));
            numsync = bitand(tmpy,65535);

            valid     = [valid;    (chan~=1)&(chan~=2)&(chan~=15)];
            overflow  = [overflow; (chan==15&tmpt==0)];
            flag      = [flag;      chan + uint32((chan==15)&(tmpt~=0))];

            mark = find(flag==16);


            cs_ovl = start_time + 2^16*cumsum(overflow);
            cs_ovl(1:end-numel(numsync)) = [];

            y      = [y; sync*(double(numsync) + cs_ovl)];
            tcspc  = [tcspc; double(tmpt)]; 

        end;

        start = mark(1);
        stop  = mark(5);

        y(1)   = [];

        if (num>0)
            if (length(y)>stop)
                y(stop+1:end)        = [];
                tcspc(stop+1:end)    = [];
                valid(stop+1:end)    = [];
                overflow(stop+1:end) = [];
                flag(stop+1:end)     = [];
            end;

            num       = stop;
            cs_ovl    = start_time + 2^16*cumsum(overflow);
            if length(cs_ovl)>0
                overcount = cs_ovl(end);
            end;

            if start>1
                y(1:start-1)        = [];
                tcspc(1:start-1)    = [];
                valid(1:start-1)    = [];
                overflow(1:start-1) = [];
                flag(1:start-1)     = [];
            end;

            y(valid==1)        = [];
            tcspc(valid==1)    = [];
            flag(valid==1)     = [];

            y(overflow==1)     = [];
            tcspc(overflow==1) = [];
            flag(overflow==1)  = [];

        end;

    end
    fclose(fin);
else
    disp('Not a pt3-file!');
end
