function [tt, ngbr] = FErode(t, c1, c2, m);

% tout = FErode(tin, c, m); erodes the binary image t m times with coefficient c
% (C) Joerg Enderlein, 2001

tt = t;
col = zeros(size(t,1),1);
row = zeros(1,size(t,2));
col1 = zeros(size(t,1)-1,1);
row1 = zeros(1,size(t,2)-1);
% erosion
for j=1:m
   ind = tt==0;
   ngbr = 3*([ind(2:end,:); row] + [row; ind(1:end-1,:)] + [ind(:,2:end) col] + [col ind(:,1:end-1)]);
   ngbr = ngbr + 2*([[ind(2:end,2:end); row1] col] + [col [row1; ind(1:end-1,1:end-1)]] + ...
      [col [ind(2:end,1:end-1); row1]] + [[row1; ind(1:end-1,2:end)] col]);
   mim(ngbr); colorbar; drawnow; pause
   tt(ngbr>=c1 & ngbr<=c2) = 0;
end   


