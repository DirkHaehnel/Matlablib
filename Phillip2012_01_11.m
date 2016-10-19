close all

dirname = 'm:\MTusers\Phillip\120111\';
fnames = dir([dirname '*PBS max laserpower2.ht3']);
%fnames = dir([dirname 'hPinAlexa488Alexa647 500pM 7m GdCl max laserpower.ht3']);
%dirname = 'm:\MTusers\Phillip\120109\';
%fnames = dir([dirname '*PBS.ht3']);

laser_b1 = 1;
laser_b2 = 2;
laser_r1 = 3;
laser_r2 = 4;
detector_r1 = 0;
detector_r2 = 4;
detector_b1 = 8;
detector_b2 = 12;
cnum = 2;

for jf=1:1
    [y, x, chan, special, num, head] = ht3v2read([dirname fnames(jf).name], [1 1e5]);
    Resolution = head.Resolution;
    dind = unique(chan);
    dnum = length(dind);
    
    % Read TCSPC Data
    [bin, tcspcdata] = ReadTCSPC([dirname fnames(jf).name]);
    
    % Detect Laser Puls Time Gates
    [t1, len] = AutodetectTimeGates(tcspcdata, cnum);
    len = len - 20;
    t2 = 1+t1-t1(1);
    bin = [bin(t1:end); bin(1:t1(1)-1)];
    tcspcdata = [tcspcdata(t1:end,:); tcspcdata(1:t1(1)-1,:)];
    
    tcspc = zeros(len,cnum,dnum);
    tau = tcspc;
    for j=1:cnum
        for k=1:dnum
            tau(1:len,j,k) = bin(t2(j):t2(j)+len-1);
            tcspc(1:len,j,k) = tcspcdata(t2(j):t2(j)+len-1,k);
        end
    end
    tcspc(any(any(tau==0,2),3),:,:)=[];
    tau(any(any(tau==0,2),3),:,:)=[];
    semilogy(Resolution*reshape(tau,[len cnum*dnum]),reshape(tcspc,[len cnum*dnum])); drawnow
    
    %     res.tau = Resolution*tau;
    %     res.tcspc = tcspc;
    %     res.bin = bin;
    %     res.tcspcdata = tcspcdata;
    
    figure
    tst = 1;
    cnt = 0;
    bin = 2.5e3;
    num = 0:120;
    h2 = zeros(numel(num));
    while tst>0
        [tmpy, tmpx, flv, bla, tst] = ht3v2read([dirname fnames(jf).name], [cnt+1, 1e6]);
        cnt = cnt + tst;
        if (tst>0) && (~isempty(tmpy))
            x = tttr2bin(tmpy, bin, flv);
            %plot([sum(x(:,1:2),2) sum(x(:,3:4),2)]);
            %plot(sum(x(:,1:2),2), sum(x(:,3:4),2),'o')
            h2 = h2 + mHist2(sum(x(:,1:2),2), sum(x(:,3:4),2), num, num);
            pcolor(num,num,log(h2)); shading flat; axis image
            drawnow;
        end
    end
end