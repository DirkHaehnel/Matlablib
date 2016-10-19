function [y, t] = FCSCrossRead(name,t,flag,filt)

if ischar(name)
    load(name);
else
    res = name;
end

if nargin>1 && ~isempty(t)
    ind = res.autotime>=min(t) & res.autotime<=max(t);
else
    ind = res.autotime>0;
end

if nargin > 2
    if ischar(flag) && (strcmp(flag, 'normalize'))
        ndet = size(res.auto, 2);
        for jj = 1:size(res.auto, 4)
            for det1 = 1:ndet
                for det2 = 1:ndet
                    res.auto(:,det1,det2,jj) = ((res.auto(:,det1,det2,jj) / (res.rate(jj,det1) * res.rate(jj,det2)))) - 1;
                end
            end
        end
    end
end


Flag_4L2D = 0;
if nargin>2 && ~isempty(flag)
    if ischar(flag) && (strcmp(flag, '4L2D'))
        Flag_4L2D = 1;
    else
    y(:,1,:) = res.auto(ind,1,1,:)+res.auto(ind,2,2,:)+res.auto(ind,1,2,:)+res.auto(ind,2,1,:);
    y(:,2,:) = res.auto(ind,3,3,:)+res.auto(ind,4,4,:)+res.auto(ind,3,4,:)+res.auto(ind,4,3,:);
    y(:,3,:) = res.auto(ind,1,4,:)+res.auto(ind,2,3,:)+res.auto(ind,1,3,:)+res.auto(ind,2,4,:);
    y(:,4,:) = res.auto(ind,4,1,:)+res.auto(ind,3,2,:)+res.auto(ind,3,1,:)+res.auto(ind,4,2,:);
    end
end

if (size(res.auto, 2) == 16)        
    % 4 Laser (Focus) -> 4 Detector
    laser_b1 = 1;
    laser_b2 = 2;
    laser_r1 = 3;
    laser_r2 = 4;
    detector_r1 = 0;
    detector_r2 = 4;
    detector_b1 = 8;
    detector_b2 = 12;

    % red ACF + CCF
    
    % ACF 1:    red->red focus red 1->1
    y(:,1,:) = res.auto(ind, laser_r1+detector_r1, laser_r1+detector_r2, :) + res.auto(ind, laser_r1+detector_r2, laser_r1+detector_r1, :);
    % ACF 2:    red->red focus red 2->2
    y(:,2,:) = res.auto(ind, laser_r2+detector_r1, laser_r2+detector_r2, :) + res.auto(ind, laser_r2+detector_r2, laser_r2+detector_r1, :);
    % CCF 12:   red->red focus red 1->2
    y(:,3,:) = res.auto(ind, laser_r1+detector_r1, laser_r2+detector_r2, :);
    y(:,5,:) = res.auto(ind, laser_r1+detector_r2, laser_r2+detector_r1, :);
    % CCF 21:   red->red focus red 2->1
    y(:,4,:) = res.auto(ind, laser_r2+detector_r1, laser_r1+detector_r2, :);
    y(:,6,:) = res.auto(ind, laser_r2+detector_r2, laser_r1+detector_r1, :);
    
    % blue ACF + CCF

    % ACF 3:    blue->blue focus blue 1->1
    y(:,7,:) = res.auto(ind, laser_b1+detector_b1, laser_b1+detector_b2, :) + res.auto(ind, laser_b1+detector_b2, laser_b1+detector_b1, :);
    % ACF 4:    blue->blue focus blue 2->2
    y(:,8,:) = res.auto(ind, laser_b2+detector_b1, laser_b2+detector_b2, :) + res.auto(ind, laser_b2+detector_b2, laser_b2+detector_b1, :);
    % CCF 34:   blue->blue focus blue 1->2
    y(:,9,:)  = res.auto(ind, laser_b1+detector_b1, laser_b2+detector_b2, :);
    y(:,11,:) = res.auto(ind, laser_b1+detector_b2, laser_b2+detector_b1, :);
    % CCF 43:   blue->blue focus blue 2->1
    y(:,10,:) = res.auto(ind, laser_b2+detector_b1, laser_b1+detector_b2, :);
    y(:,12,:) = res.auto(ind, laser_b2+detector_b2, laser_b1+detector_b1, :);

    % FRET ACF + CCF

    % ACF 5:    red->red focus blue 1->1
    y(:,13,:) = res.auto(ind, laser_b1+detector_r1, laser_b1+detector_r2, :) + res.auto(ind, laser_b1+detector_r2, laser_b1+detector_r1, :);
    % ACF 6:    red->red focus blue 2->2
    y(:,14,:) = res.auto(ind, laser_b2+detector_r1, laser_b2+detector_r2, :) + res.auto(ind, laser_b2+detector_r2, laser_b2+detector_r1, :);
    % CCF 56:   red->red focus blue 1->2
    y(:,15,:) = res.auto(ind, laser_b1+detector_r1, laser_b2+detector_r2, :);
    y(:,17,:) = res.auto(ind, laser_b1+detector_r2, laser_b2+detector_r1, :);
    % CCF 65:   red->red focus blue 2->1
    y(:,16,:) = res.auto(ind, laser_b2+detector_r1, laser_b1+detector_r2, :);
    y(:,18,:) = res.auto(ind, laser_b2+detector_r2, laser_b1+detector_r1, :);

    % 2 Color CCF + 2 Color - 2 Focus CCCF

    % ACF 7:    red<->blue focus 1->1
    tmp(:,1,:) = res.auto(ind, laser_b1+detector_b1, laser_r1+detector_r2, :) + res.auto(ind, laser_r1+detector_r2, laser_b1+detector_b1, :);
    tmp(:,2,:) = res.auto(ind, laser_b1+detector_b2, laser_r1+detector_r1, :) + res.auto(ind, laser_r1+detector_r1, laser_b1+detector_b2, :);
    tmp(:,3,:) = res.auto(ind, laser_b1+detector_b1, laser_r1+detector_r1, :) + res.auto(ind, laser_r1+detector_r1, laser_b1+detector_b1, :);
    tmp(:,4,:) = res.auto(ind, laser_b1+detector_b2, laser_r1+detector_r2, :) + res.auto(ind, laser_r1+detector_r2, laser_b1+detector_b2, :);
    y(:,19,:) = (tmp(:,1,:) + tmp(:,2,:) + tmp(:,3,:) + tmp(:,4,:)) / 4;
    % ACF 8:    red<->blue focus 2->2
    tmp(:,1,:) = res.auto(ind, laser_b2+detector_b1, laser_r2+detector_r2, :) + res.auto(ind, laser_r2+detector_r2, laser_b2+detector_b1, :);
    tmp(:,2,:) = res.auto(ind, laser_b2+detector_b2, laser_r2+detector_r1, :) + res.auto(ind, laser_r2+detector_r1, laser_b2+detector_b2, :);
    tmp(:,3,:) = res.auto(ind, laser_b2+detector_b1, laser_r2+detector_r1, :) + res.auto(ind, laser_r2+detector_r1, laser_b2+detector_b1, :);
    tmp(:,4,:) = res.auto(ind, laser_b2+detector_b2, laser_r2+detector_r2, :) + res.auto(ind, laser_r2+detector_r2, laser_b2+detector_b2, :);
    y(:,20,:) = (tmp(:,1,:) + tmp(:,2,:) + tmp(:,3,:) + tmp(:,4,:)) / 4;
    % CCF 78:   red<->blue focus 1->2
    tmp(:,1,:) = res.auto(ind, laser_b1+detector_b1, laser_r2+detector_r2, :) + res.auto(ind, laser_r1+detector_r2, laser_b2+detector_b1, :);
    tmp(:,2,:) = res.auto(ind, laser_b1+detector_b2, laser_r2+detector_r1, :) + res.auto(ind, laser_r1+detector_r1, laser_b2+detector_b2, :);
    tmp(:,3,:) = res.auto(ind, laser_b1+detector_b1, laser_r2+detector_r1, :) + res.auto(ind, laser_r1+detector_r1, laser_b2+detector_b1, :);
    tmp(:,4,:) = res.auto(ind, laser_b1+detector_b2, laser_r2+detector_r2, :) + res.auto(ind, laser_r1+detector_r2, laser_b2+detector_b2, :);
    y(:,21,:) = (tmp(:,1,:) + tmp(:,2,:) + tmp(:,3,:) + tmp(:,4,:)) / 4;
    % CCF 87:   red<->blue focus 2->1
    tmp(:,1,:) = res.auto(ind, laser_b2+detector_b1, laser_r1+detector_r2, :) + res.auto(ind, laser_r2+detector_r2, laser_b1+detector_b1, :);
    tmp(:,2,:) = res.auto(ind, laser_b2+detector_b2, laser_r1+detector_r1, :) + res.auto(ind, laser_r2+detector_r1, laser_b1+detector_b2, :);
    tmp(:,3,:) = res.auto(ind, laser_b2+detector_b1, laser_r1+detector_r1, :) + res.auto(ind, laser_r2+detector_r1, laser_b1+detector_b1, :);
    tmp(:,4,:) = res.auto(ind, laser_b2+detector_b2, laser_r1+detector_r2, :) + res.auto(ind, laser_r2+detector_r2, laser_b1+detector_b2, :);   
    y(:,22,:) = (tmp(:,1,:) + tmp(:,2,:) + tmp(:,3,:) + tmp(:,4,:)) / 4;

    % donor vs acceptor cross-correlation upon donor excitation

    % ACF:    focus 1->1
    y(:,23,:) = res.auto(ind, laser_b1+detector_b1, laser_b1+detector_r1, :) + res.auto(ind, laser_b1+detector_b1, laser_b1+detector_r2, :) + ...
        res.auto(ind, laser_b1+detector_b2, laser_b1+detector_r1, :) + res.auto(ind, laser_b1+detector_b2, laser_b1+detector_r2, :);
    % ACF:    focus 2->2
    y(:,24,:) = res.auto(ind, laser_b2+detector_b1, laser_b2+detector_r1, :) + res.auto(ind, laser_b2+detector_b1, laser_b2+detector_r2, :) + ...
        res.auto(ind, laser_b2+detector_b2, laser_b2+detector_r1, :) + res.auto(ind, laser_b2+detector_b2, laser_b2+detector_r2, :);
    % CCF:   focus 1->2
    y(:,25,:) = (res.auto(ind, laser_b1+detector_b1, laser_b2+detector_r1, :) + res.auto(ind, laser_b1+detector_b1, laser_b2+detector_r2, :) + ...
        res.auto(ind, laser_b1+detector_b2, laser_b2+detector_r1, :) + res.auto(ind, laser_b1+detector_b2, laser_b2+detector_r2, :));
    % CCF:   focus 2->1
    y(:,26,:) = (res.auto(ind, laser_b2+detector_b1, laser_b1+detector_r1, :) + res.auto(ind, laser_b2+detector_b1, laser_b1+detector_r2, :) + ...
        res.auto(ind, laser_b2+detector_b2, laser_b1+detector_r1, :) + res.auto(ind, laser_b2+detector_b2, laser_b1+detector_r2, :));
    
elseif (size(res.auto, 2) == 8)
    if (Flag_4L2D == 1)
        % 4 Laser (Focus) -> 2 Detector
        laser_b1 = 1;
        laser_b2 = 2;
        laser_r1 = 3;
        laser_r2 = 4;
        detector_r1 = 0;
        detector_r2 = 4;

        % RED ACF + CCF

        % ACF 1:    blue->blue focus blue 1->1
        y(:,1,:) = res.auto(ind, laser_r1+detector_r1, laser_r1+detector_r2, :) + res.auto(ind, laser_r1+detector_r2, laser_r1+detector_r1, :);
        % ACF 2:    blue->blue focus blue 2->2
        y(:,2,:) = res.auto(ind, laser_r2+detector_r1, laser_r2+detector_r2, :) + res.auto(ind, laser_r2+detector_r2, laser_r2+detector_r1, :);
        % CCF 12:   blue->blue focus blue 1->2
        y(:,3,:) = res.auto(ind, laser_r1+detector_r1, laser_r2+detector_r2, :);
        y(:,5,:) = res.auto(ind, laser_r1+detector_r2, laser_r2+detector_r1, :);
        % CCF 21:   blue->blue focus blue 2->1
        y(:,4,:) = res.auto(ind, laser_r2+detector_r1, laser_r1+detector_r2, :);
        y(:,6,:) = res.auto(ind, laser_r2+detector_r2, laser_r1+detector_r1, :);

        % FRET ACF + CCF

        % ACF 3:    red->red focus blue 1->1
        y(:,7,:) = res.auto(ind, laser_b1+detector_r1, laser_b1+detector_r2, :) + res.auto(ind, laser_b1+detector_r2, laser_b1+detector_r1, :);
        % ACF 4:    red->red focus blue 2->2
        y(:,8,:) = res.auto(ind, laser_b2+detector_r1, laser_b2+detector_r2, :) + res.auto(ind, laser_b2+detector_r2, laser_b2+detector_r1, :);
        % CCF 34:   red->red focus blue 1->2
        y(:,9,:)  = res.auto(ind, laser_b1+detector_r1, laser_b2+detector_r2, :);
        y(:,11,:) = res.auto(ind, laser_b1+detector_r2, laser_b2+detector_r1, :);
        % CCF 43:   red->red focus blue 2->1
        y(:,10,:) = res.auto(ind, laser_b2+detector_r1, laser_b1+detector_r2, :);
        y(:,12,:) = res.auto(ind, laser_b2+detector_r2, laser_b1+detector_r1, :);

    else
        % 2 Laser (Focus) -> 4 Detector
        laser_b1 = 1;
        laser_b2 = 2;
        detector_r1 = 0;
        detector_r2 = 2;
        detector_b1 = 4;
        detector_b2 = 6;

        % blue ACF + CCF

        % ACF 1:    blue->blue focus blue 1->1
        y(:,1,:) = res.auto(ind, laser_b1+detector_b1, laser_b1+detector_b2, :) + res.auto(ind, laser_b1+detector_b2, laser_b1+detector_b1, :);
        % ACF 2:    blue->blue focus blue 2->2
        y(:,2,:) = res.auto(ind, laser_b2+detector_b1, laser_b2+detector_b2, :) + res.auto(ind, laser_b2+detector_b2, laser_b2+detector_b1, :);
        % CCF 12:   blue->blue focus blue 1->2
        y(:,3,:) = res.auto(ind, laser_b1+detector_b1, laser_b2+detector_b2, :);
        y(:,5,:) = res.auto(ind, laser_b1+detector_b2, laser_b2+detector_b1, :);
        % CCF 21:   blue->blue focus blue 2->1
        y(:,4,:) = res.auto(ind, laser_b2+detector_b1, laser_b1+detector_b2, :);
        y(:,6,:) = res.auto(ind, laser_b2+detector_b2, laser_b1+detector_b1, :);

        % FRET ACF + CCF

        % ACF 3:    red->red focus blue 1->1
        y(:,7,:) = res.auto(ind, laser_b1+detector_r1, laser_b1+detector_r2, :) + res.auto(ind, laser_b1+detector_r2, laser_b1+detector_r1, :);
        % ACF 4:    red->red focus blue 2->2
        y(:,8,:) = res.auto(ind, laser_b2+detector_r1, laser_b2+detector_r2, :) + res.auto(ind, laser_b2+detector_r2, laser_b2+detector_r1, :);
        % CCF 34:   red->red focus blue 1->2
        y(:,9,:)  = res.auto(ind, laser_b1+detector_r1, laser_b2+detector_r2, :);
        y(:,11,:) = res.auto(ind, laser_b1+detector_r2, laser_b2+detector_r1, :);
        % CCF 43:   red->red focus blue 2->1
        y(:,10,:) = res.auto(ind, laser_b2+detector_r1, laser_b1+detector_r2, :);
        y(:,12,:) = res.auto(ind, laser_b2+detector_r2, laser_b1+detector_r1, :);

        % 2 Color CCF + 2 Color - 2 Focus CCCF

        % ACF 5:    red<->blue focus 1->1
        tmp(:,1,:) = res.auto(ind, laser_b1+detector_b1, laser_b1+detector_r2, :) + res.auto(ind, laser_b1+detector_r2, laser_b1+detector_b1, :);
        tmp(:,2,:) = res.auto(ind, laser_b1+detector_b2, laser_b1+detector_r1, :) + res.auto(ind, laser_b1+detector_r1, laser_b1+detector_b2, :);
        tmp(:,3,:) = res.auto(ind, laser_b1+detector_b1, laser_b1+detector_r1, :) + res.auto(ind, laser_b1+detector_r1, laser_b1+detector_b1, :);
        tmp(:,4,:) = res.auto(ind, laser_b1+detector_b2, laser_b1+detector_r2, :) + res.auto(ind, laser_b1+detector_r2, laser_b1+detector_b2, :);
        y(:,13,:) = (tmp(:,1,:) + tmp(:,2,:) + tmp(:,3,:) + tmp(:,4,:)) / 4;
        % ACF 6:    red<->blue focus 2->2
        tmp(:,1,:) = res.auto(ind, laser_b2+detector_b1, laser_b2+detector_r2, :) + res.auto(ind, laser_b2+detector_r2, laser_b2+detector_b1, :);
        tmp(:,2,:) = res.auto(ind, laser_b2+detector_b2, laser_b2+detector_r1, :) + res.auto(ind, laser_b2+detector_r1, laser_b2+detector_b2, :);
        tmp(:,3,:) = res.auto(ind, laser_b2+detector_b1, laser_b2+detector_r1, :) + res.auto(ind, laser_b2+detector_r1, laser_b2+detector_b1, :);
        tmp(:,4,:) = res.auto(ind, laser_b2+detector_b2, laser_b2+detector_r2, :) + res.auto(ind, laser_b2+detector_r2, laser_b2+detector_b2, :);
        y(:,14,:) = (tmp(:,1,:) + tmp(:,2,:) + tmp(:,3,:) + tmp(:,4,:)) / 4;
        % CCF 56:   red<->blue focus 1->2
        tmp(:,1,:) = res.auto(ind, laser_b1+detector_b1, laser_b2+detector_r2, :) + res.auto(ind, laser_b1+detector_r2, laser_b2+detector_b1, :);
        tmp(:,2,:) = res.auto(ind, laser_b1+detector_b2, laser_b2+detector_r1, :) + res.auto(ind, laser_b1+detector_r1, laser_b2+detector_b2, :);
        tmp(:,3,:) = res.auto(ind, laser_b1+detector_b1, laser_b2+detector_r1, :) + res.auto(ind, laser_b1+detector_r1, laser_b2+detector_b1, :);
        tmp(:,4,:) = res.auto(ind, laser_b1+detector_b2, laser_b2+detector_r2, :) + res.auto(ind, laser_b1+detector_r2, laser_b2+detector_b2, :);
        y(:,15,:) = (tmp(:,1,:) + tmp(:,2,:) + tmp(:,3,:) + tmp(:,4,:)) / 4;
        % CCF 65:   red<->blue focus 2->1
        tmp(:,1,:) = res.auto(ind, laser_b2+detector_b1, laser_b1+detector_r2, :) + res.auto(ind, laser_b2+detector_r2, laser_b1+detector_b1, :);
        tmp(:,2,:) = res.auto(ind, laser_b2+detector_b2, laser_b1+detector_r1, :) + res.auto(ind, laser_b2+detector_r1, laser_b1+detector_b2, :);
        tmp(:,3,:) = res.auto(ind, laser_b2+detector_b1, laser_b1+detector_r1, :) + res.auto(ind, laser_b2+detector_r1, laser_b1+detector_b1, :);
        tmp(:,4,:) = res.auto(ind, laser_b2+detector_b2, laser_b1+detector_r2, :) + res.auto(ind, laser_b2+detector_r2, laser_b1+detector_b2, :);   
        y(:,16,:) = (tmp(:,1,:) + tmp(:,2,:) + tmp(:,3,:) + tmp(:,4,:)) / 4;
    end
    
else
    % 2 Laser -> 2 Detector
    L1 = 1;
    L2 = 2;
    D1 = 0;
    D2 = 2;

    y(:,1,:) = res.auto(ind,L1+D2,L1+D2,:);
    y(:,2,:) = res.auto(ind,L2+D1,L2+D1,:);
    y(:,3,:) = res.auto(ind,L1+D2,L2+D1,:) / 2;
    y(:,4,:) = res.auto(ind,L2+D1,L1+D2,:) / 2;
    y(:,5,:) = res.auto(ind,L1+D2,L2+D1,:) / 2;
    y(:,6,:) = res.auto(ind,L2+D1,L1+D2,:) / 2;
    
    y(:,7,:) = res.auto(ind,L1+D1,L1+D1,:);
    y(:,8,:) = res.auto(ind,L2+D2,L2+D2,:);
    y(:,9,:) = res.auto(ind,L1+D1,L2+D2,:) / 2;
    y(:,10,:) = res.auto(ind,L2+D2,L1+D1,:) / 2;
    y(:,11,:) = res.auto(ind,L1+D1,L2+D2,:) / 2;
    y(:,12,:) = res.auto(ind,L2+D2,L1+D1,:) / 2;    
end

if nargin>3 && ~isempty(filt)
    p = Simplex('ExpFun',3,[],[],[],[],1:size(res.rate,1),sum(res.rate,2));
    [bla, bla, bla, z] = ExpFun(p,1:size(res.rate,1),sum(res.rate,2));
    y = y(:,:,sum(res.rate,2)<z);
end

t = res.autotime(ind);
l = length(res.time);
threshold = mean(res.time(1:l-1));
if (res.time(l) < 0.6 * threshold)
    y(:,:,l) = [];
end

