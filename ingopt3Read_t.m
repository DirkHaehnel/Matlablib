function [y, tcspc, flag, num, overcount, head] = pt3Read_t(name, start_cnt, start_time, max_time, sync)

%        [y, tcspc, flag, num, overcount, head] = pt3Read_t(name, start_cnt, start_time, max_time, sync)


if strcmp(name(end-3:end),'.pt3')
    hlength = 0;

    if nargin>1

        if (nargin<5)||(nargout==6) 
            [head, hlength] = pt3Read_head(name);
            sync            = 1e9/head.CntRate0;
        end;

        if nargin<4
            max_time = 0;
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
        %         if (hlength==0)
        %             hlength = fread(fin, 1, 'int32');
        %             if hlength > 0
        %                 fclose(fin);
        %                 pt3_headerupdate(name);
        %                 fin = fopen(name,'r');
        %                 hlength = fread(fin, 1, 'int32');
        %             end;
        %             hlength = abs(hlength);
        %         end

        fseek(fin, hlength+4*(start_cnt-1), 'bof');

        valid     = [];
        overflow  = [];
        y         = 0;
        tcspc     = [];
        flag      = [];

        num       = 0;
        overcount = start_time;

        while (y(end)<max_time)&&(~feof(fin))

            [tmpy num] = fread(fin, 5E5, 'uint32');

            tmpy = uint32(tmpy);

            chan    = bitshift(tmpy, -28);
            tmpt    = double(bitshift(bitand(tmpy,268369920),-16));
            numsync = bitand(tmpy,65535);
            
            valid     = [valid;    (chan~=1)&(chan~=2)&(chan~=3)&(chan~=15)];
            overflow  = [overflow; (chan==15&tmpt==0)];
            flag      = [flag;      chan + uint32((chan==15)&(tmpt~=0))];

            cs_ovl                = start_time + 2^16*cumsum(overflow);
            cs_ovl(1:length(y)-1) = [];

            y     = [y; sync*(double(numsync) + cs_ovl)];
            tcspc = [tcspc; double(tmpt)]; 

        end;
        
        y(1)   = [];
        ind    = (y>=max_time);
        
        if (num>0)
            y(ind)        = [];
            tcspc(ind)    = [];
            valid(ind)    = [];
            overflow(ind) = [];
            flag(ind)     = [];

            num       = length(y);
            cs_ovl    = start_time + 2^16*cumsum(overflow);
            if length(cs_ovl)>0
                overcount = cs_ovl(end);
            end;

            y(valid==1)     = [];
            tcspc(valid==1) = [];
            flag(valid==1)  = [];
        end;

    end
    fclose(fin);
    
else
    disp('Not a pt3-file!');
end
