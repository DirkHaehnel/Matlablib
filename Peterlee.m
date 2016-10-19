% Program Peterlee
% PETERLEE('Name') reads SMD data, applies a Lee filter,
% and sifts for bursts
%
% bs gives the burst sizes
% bh gives the burst heights
% trans gives the transit times

fid=fopen('ma971217.01')
tmp = setstr(fread(fid, 512, 'char'))';
N = str2num(tmp(findstr(tmp, '$TOT')+5:findstr(tmp, '$N_SCAN')-2));
scans = str2num(tmp(findstr(tmp, '$N_SCAN')+8:findstr(tmp, '$DWELL')-2))

if N ==[] 
	N = 2^13;
else
	N = N/scans
end

lee = 5; % width of Lee filter
maxheight = 30; % maximum expected height of counts per bin
bs = [];
bh = [];
trans = [];
bsd = [];
bhd = [];
ttd = [];
t = 1:N;
bstmp = zeros(1,N);
bhtmp = zeros(1,N);
mn = zeros(1,N);
sig = zeros(1,N);
tmp = zeros(2*lee+1, N-2*lee);
sig0 = 5;
close all
for j = 1:12
	j
	x = fread(fid, N, 'int32')';
end
for j = 1:100
	j
	x = fread(fid, N, 'int32')';
	bg = sum(x==1)/sum(x==0)
	tsh = 1.5*bg;
	for k = 1:lee
		mn(k) = mean(x(1:k+lee));
	end
	for k=1:2*lee+1
		tmp(k,:) = x(k:N-2*lee-1+k);
	end;
	mn(lee+1:N-lee) = mean(tmp);
	for k = N-lee+1:N
		mn(k) = mean(x(k-lee:N));
	end
	for k = 1:lee
	sig(k) = mean((x(1:k+lee)-mn(1:k+lee)).^2);
	end
	for k=1:2*lee+1
		tmp(k,:) = (x(k:N-2*lee-1+k) - mn(k:N-2*lee-1+k)).^2;
	end;
	sig(lee+1:N-lee) = mean(tmp);
	for k = N-lee+1:N
	sig(k) = mean((x(k-lee:N)-mn(k-lee:N)).^2);
	end
	y = mn + (x-mn).*sig./(sig + sig0);
	plot(t,y,[1 N],[tsh tsh]); drawnow;
	t1=t([0 y(1:N-1)]<=tsh & y>tsh);
	t2=t([y(2:N) 0]<=tsh & y>tsh);
	for k = 1:length(t1)
		bstmp(k) = sum(x(t1(k):t2(k))) - (t2(k)-t1(k)+1)*bg;
		bhtmp(k) = max(x(t1(k):t2(k)));
	end
	bs = [bs bstmp(1:length(t1))];
%	plot(bs); drawnow;
	bh = [bh bhtmp(1:length(t1))];
	trans = [trans (t2-t1)+1];
	if rem(j,10)==0
		if bsd==[]
			mbs = max(bs);
			mbh = max(bh);
			mtr = max(trans);
		bsd = hist(bs, 1:mbs);
			bhd = hist(bh, 1:mbh);
			ttd = hist(trans, 1:mtr);
		else
			bsd = [bsd; hist(bs, 1:mbs)];
			bhd = [bhd; hist(bh, 1:mbh)];
			ttd = [ttd; hist(trans, 1:mtr)];
		end
		bs = [];
		bh = [];
		trans = [];
	end
end
fclose(fid);
