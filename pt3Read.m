function [y, tcspc, flag, num, overcount, head] = pt3Read(name, start_cnt, cnts, start_time, sync)

%        [y, tcspc, flag, num, overcount, head] = pt3Read(name, start_cnt, cnts, start_time, sync);

if strcmp(name(end-3:end),'.pt3')
    
    hlength = 0;
    
    if nargin>1
        if (nargin<5)||(nargout==6) 
            [head, hlength] = pt3Read_head(name);
            sync      = 1e9/head.CntRate0;
        end;
        if nargin<4
            start_time = 0;
        end;

        if nargin<3
            cnts = 0;
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

        overcount = start_time;

        [tmpy num] = fread(fin, cnts, 'uint32');
        
        tmpy = uint32(tmpy);

        chan    = bitshift(tmpy, -28);
        tcspc   = double(bitshift(bitand(tmpy,268369920),-16));
        numsync = bitand(tmpy,65535);

        valid     = (chan<1)|(chan>4)&(chan<15);
        overflow  = (chan==15)&(tcspc==0);
        flag      = chan + uint32((chan==15)&(tcspc~=0));
        
        cs_ovl    = start_time + 2^16*cumsum(overflow);
        y         = sync*(double(numsync) + cs_ovl);

        if (num>0)

            if ~isempty(cs_ovl)
                overcount = cs_ovl(end);
            end;

            y(valid==1 | overflow==1)        = [];
            tcspc(valid==1 | overflow==1)    = [];
            flag(valid==1 | overflow==1)     = [];

        end;

    end
    fclose(fin);

else
    disp('Not a pt3-file!');
end


