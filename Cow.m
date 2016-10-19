function [y, cx, cy] = cow(x,m);

% function [y, cx, cy] = cow(x) calculates the center of weight image y for image x with gyration radius m

mm = disk(m); mm = mm/sum(sum(mm));
[jx,jy] = meshgrid(-m+0.5:size(x,1)+m,-m+0.5:size(x,2)+m);
jx = jx'; jy = jy';
tmp = [x(m:-1:1,m:-1:1) x(m:-1:1,:) x(m:-1:1,end:-1:end-m+1); ...
      x(:,m:-1:1) x x(:,end:-1:end-m+1); ...
      x(end:-1:end-m+1,m:-1:1) x(end:-1:end-m+1,:) x(end:-1:end-m+1,end:-1:end-m+1)];
jx = jx.*tmp;
jy = jy.*tmp;
tmp = mconv2(x,mm);
cx = mconv2(jx,mm);
cx = cx(m+1:end-m,m+1:end-m)./tmp;
cy = mconv2(jy,mm);
cy = cy(m+1:end-m,m+1:end-m)./tmp;
y = mhist3(reshape(cx,prod(size(cx)),1),reshape(cy,prod(size(cy)),1),0.5:size(cx,1),0.5:size(cx,2))';