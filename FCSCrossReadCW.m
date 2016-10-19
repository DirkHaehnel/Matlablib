function [y, t] = FCSCrossReadCW(name)

if ischar(name)
    load(name);
else
    res = name;
end

if (size(res.auto, 2) == 2)
    % 2 Detector
    detector_p1 = 1;
    detector_o1 = 2;
    % CCF 12:  parallel -> orthogonal
    y(:,1,:) = res.auto(:, detector_p1, detector_o1, :);
    % CCF 21:  orthogonal -> parallel
    y(:,2,:) = res.auto(:, detector_o1, detector_p1, :);
elseif (size(res.auto, 2) == 4)        
    % 4 Detector
    detector_p1 = 1;
    detector_p2 = 2;
    detector_o1 = 3;
    detector_o2 = 4;

    % red ACF + CCF
    
    % ACF 1:   parallel -> parallel
    y(:,1,:) = (res.auto(:, detector_p1, detector_p2, :) + res.auto(:, detector_p2, detector_p1, :)) / 2;
    % ACF 2:   orthogonal -> orthogonal
    y(:,2,:) = (res.auto(:, detector_o1, detector_o2, :) + res.auto(:, detector_o2, detector_o1, :)) / 2;
    % CCF 12:  parallel -> orthogonal
    tmp(:,1,:) = res.auto(:, detector_p1, detector_o1, :) + res.auto(:, detector_p1, detector_o2, :);
    tmp(:,2,:) = res.auto(:, detector_p2, detector_o1, :) + res.auto(:, detector_p2, detector_o2, :);
    y(:,3,:) = (tmp(:,1,:) + tmp(:,2,:)) / 4;
    % CCF 21:  orthogonal -> parallel
    tmp(:,1,:) = res.auto(:, detector_o1, detector_p1, :) + res.auto(:, detector_o1, detector_p2, :);
    tmp(:,2,:) = res.auto(:, detector_o2, detector_p1, :) + res.auto(:, detector_o2, detector_p2, :);
    y(:,4,:) = (tmp(:,1,:) + tmp(:,2,:)) / 4;
else
    y = [];
end

t = res.autotime;
