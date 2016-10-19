function [res, head] = PIE2FRET(name)

photons = 1e7;

close all

if strcmp(name(end-2:end),'pt3')
    head = pt3v2Read(name);
    Resolution = head.Resolution;
    NChannels  = 2^12;
    dind = [1 3];
elseif strcmp(name(end-2:end),'t3r')
    head = tttrRead(name);
    Resolution = head.Resolution;
    NChannels  = 2^12;
	dind = [0 1];
elseif strcmp(name(end-2:end),'ht3')
    if nargin<3 || isempty(para)
        Timeunit = 50e-9;
        Resolution = 2e-3;
    else
        Timeunit = para(1)*1e-9;
        Resolution = para(2);
    end
    NChannels  = ceil(Timeunit/(1e-9*Resolution));
    dind = [0 1];
else
    disp('Not a valid file!');
    return
end

tcspcdata = zeros(NChannels,2);
tst = 1;
cnt = 0;
bin = 0:NChannels-1;
while tst>0
    if strcmp(name(end-2:end),'pt3')
        [tmpy, tmpx, flv, bla, tst] = pt3v2read(name, [cnt+1, photons]);
    elseif strcmp(name(end-2:end),'t3r')
        [tmpy, tmpx, flv, bla, tst] = tttrread(name, [cnt+1, photons]);
    elseif strcmp(name(end-2:end),'ht3')
        [tmpy, tmpx, flv, special, tst, overcount, head] = ht3read(name, [cnt+1, photons]);
    end
    cnt = cnt+tst;
    disp(cnt);
    if tst>0
        tcspcdata(:,1) = tcspcdata(:,1) + mHist(tmpx(flv==dind(1)),bin);
        tcspcdata(:,2) = tcspcdata(:,2) + mHist(tmpx(flv==dind(2)),bin);
    end
end

ind = sum(tcspcdata,2)==0;
bin(ind) = [];
tcspcdata(ind,:) = [];
bin([1:30 end-10:end])=[];
tcspcdata([1:30 end-10:end],:) = [];

bin = bin';
% tmp = mean(tcspcdata,2);
% [ind, num] = mCluster(tmp-min(tmp)>0.1*mean(tmp));
% t1 = max([1, min(bin(ind==length(num)))-50]);
% if length(num)>1 && num(end-1)/num(end)>0.3
%     t2 = min(bin(ind==length(num)-1))-50;
%     if t1>t2
%         tmp = t1; t1 = t2; t2 = tmp;
%     end
%     len = min([t2-t1,length(bin)-t2-20]);
%     tau(1:len,[1 3]) = bin(t1:t1+len-1)*[1 1];
%     tau(1:len,[2 4]) = bin(t2:t2+len-1)*[1 1];
%     tcspc(1:len,1) = tcspcdata(t1:t1+len-1, 1);
%     tcspc(1:len,3) = tcspcdata(t1:t1+len-1, 2);
%     tcspc(1:len,2) = tcspcdata(t2:t2+len-1, 1);
%     tcspc(1:len,4) = tcspcdata(t2:t2+len-1, 2);
%     if nargin>3 && ~isempty(cutoff)
%         ind = tau(:,1) - tau(1,1) > cutoff/Resolution;
%         tau = tau(ind,:);
%         tcspc = tcspc(ind,:);
%     end
% else
%     t2 = length(bin)-50;
%     len = t2-t1+1;
%     tau(1:len,[1 2]) = bin(t1:t2)*[1 1];
%     tcspc(1:len,1) = tcspcdata(t1:t2, 1);
%     tcspc(1:len,2) = tcspcdata(t1:t2, 2);
% end
len = floor(length(bin)/2);
tau(1:len,[1 3]) = bin(1:len)*[1 1];
tau(1:len,[2 4]) = bin(len+1:2*len)*[1 1];
tcspc(1:len,1) = tcspcdata(1:len,1);
tcspc(1:len,3) = tcspcdata(1:len,2);
tcspc(1:len,2) = tcspcdata(len+1:2*len,1);
tcspc(1:len,4) = tcspcdata(len+1:2*len,2);

tcspc(any(tau==0,2),:)=[];
tau(any(tau==0,2),:)=[];
res.tau = Resolution*tau;
res.tcspc = tcspc;
res.bin = bin;
res.tcspcdata = tcspcdata;

len = size(res.tcspc,2);
for j=1:size(res.tcspc,2)
    subplot(len/2,2,j)
    [ind ind] = max(res.tcspc(:,j));
    res.ptau{j} = res.tau(min([end-10 ind+50]):end,j);
    res.pt(j) = Simplex('ExpFun',2,0,[],[],[],res.ptau{j},res.tcspc(min([end-10 ind+50]):end,j),1);
    [bla bla bla res.z{j}] = ExpFun(res.pt(j),res.ptau{j},res.tcspc(min([end-10 ind+50]):end,j));
    drawnow; 
end
% if length(pt)>2
% 	pt(1)=Simplex('FretFun',2,0,[],[],[],res.tau(min([end-10 ind+100]):end,1),res.tcspc(min([end-10 ind+100]):end,1),[pt(3) pt(2)],1);
%     z(:,1) = FretFun(pt(1),res.tau(min([end-10 ind+100]):end,1),res.tcspc(min([end-10 ind+100]):end,1),[pt(3) pt(2)]);
% end

clf
semilogy(res.tau,res.tcspc,'o')
hold on
col=get(gca,'ColorOrder');
for j=1:length(res.ptau)
    semilogy(res.ptau{j},res.z{j},'color','y');
end
hold off
ax = axis;
if length(res.pt)>2
    text(res.tau(1,1),0.6*ax(4),mnum2str(res.pt(1),2,2),'Color',col(1,:))
    text(res.tau(1,1),0.4*ax(4),mnum2str(res.pt(3),2,2),'Color',col(3,:))    
    text(res.tau(1,2),0.6*ax(4),mnum2str(res.pt(2),2,2),'Color',col(2,:))
    text(res.tau(1,2),0.4*ax(4),mnum2str(res.pt(4),2,2),'Color',col(4,:))    
else
    text(res.tau(1,1),0.6*ax(4),mnum2str(res.pt(1),2,2),'Color',col(1,:))
    text(res.tau(1,1),0.4*ax(4),mnum2str(res.pt(2),2,2),'Color',col(2,:))
end
xlabel('time (ns)')
ylabel('counts');
