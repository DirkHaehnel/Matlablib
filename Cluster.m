function [tt, len] = Cluster(t,th);

% [tt, len] = Cluster(t,th) determines, sorts and numerates 
% clusters in an array t, whereby clusters of size smaller 
% than th are ignored. Default value th=1.
% (C) Joerg Enderlein, 2001
%
% Example:
% t=round(rand(50)-0.2); 
% subplot(121); imagesc(t); axis image; 
% subplot(122); imagesc(cluster(t,5)); axis image; 

if nargin<2
   th = 1;
end

[m,n] = size(t);
ind = (t>0);
indr = logical([zeros(m,1) ind(:,1:end-1)]);
indd = logical([zeros(1,n); ind(1:end-1,:)]);
indl = logical([ind(:,2:end) zeros(m,1)]);
indu = logical([ind(2:end,:); zeros(1,n)]);
tmpold = 0;
t = reshape(cumsum(reshape(t,1,m*n)),m,n).*t;
tmp = sum(sum(t));
while ~(tmpold==tmp)
   tmpold = tmp;
   t(indd) = min(t(indd),t(indd([2:end 1],:)));   
   t(indr) = min(t(indr),t(indr(:,[2:end 1])));
   t(indu) = min(t(indu),t(indu([end 1:end-1],:)));   
   t(indl) = min(t(indl),t(indl(:,[end 1:end-1])));
   tmp = sum(sum(t));
end
t = reshape(t,1,m*n);
[b,j,j] = unique(t);
len = hist(t(t>0),b(2:end));
[len,k] = sort(len);
tmp = sum(len>=th);
b = [zeros(1,length(b)-tmp) 1:tmp];
b([1 k+1]) = b;
tt = reshape(b(j),m,n);
len = len(len>=th);

