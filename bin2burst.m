function [bs, t1, t2] = bin2burst(x,minbs,m,cap,mm,fac)

if nargin<2 || isempty(minbs) % minimum bin size
    minbs = 10;
end
if nargin<3 || isempty(m) % half width of moving average window
    m = 10;
end
if nargin<4 || isempty(cap) % time window capping
    cap = 10;
end
if nargin<5 || isempty(mm)
    mm = 10*m;
end
if nargin<6 || isempty(fac)
    fac = 1.5;
end
wndw = 15;

xx = sum(x,2);
mn = mConv(xx, ones(2*m+1,1)/(2*m+1));
yy = mConv((xx - mn).^2, ones(2*m+1,1)/(2*m+1));

if mm==0
	bg = ones(length(xx),1);
else
	bg = mConv(xx, ones(1,2*mm+1)/(2*mm+1));
end

t = 1:length(yy);
t1=t(yy(1:end-1)<=fac*bg(1:end-1) & yy(2:end)>fac*bg(2:end)) + 1 - cap;
t2=t(yy(2:end)<=fac*bg(2:end) & yy(1:end-1)>fac*bg(1:end-1)) + cap;
t1(t1<1) = 1;
t2(t2>max(t)) = max(t);

% % background subtraction
% if mm==0
% 	bg = zeros(size(x));
% else
% 	bg = mConv(x, ones(1,2*mm+1)/(2*mm+1));
% end
% x = x - bg;

if ~(isempty(t1) || isempty(t2))
    if t2(1)<t1(1)
        t2(1)=[];
    end
    if t2(end)<t1(end)
        t1(end)=[];
    end
    if length(t2)==length(t1)+1
        t2(1)=[];
    end
    if length(t1)==length(t2)+1
        t1(end)=[];
    end

    t2(t1<wndw-1) = [];
    t1(t1<wndw-1) = [];
    t1(t2>length(xx)-wndw+1) = [];
    t2(t2>length(xx)-wndw+1) = [];

    x = cumsum(x);
    bs = x(t2,:) - x(t1-1,:);
    t1(sum(bs,2)<minbs) = [];
    t2(sum(bs,2)<minbs) = [];
    bs(sum(bs,2)<minbs,:) = [];
else
    t1 = []; t2 = []; bs = [];
end

