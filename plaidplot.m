function plaidplot(res, ind);

if nargin==1
    dt = diff(res.autotime); dt(end+1) = dt(end);
    fill3(res.autotime*0,res.nums,res.bintau,squeeze(sum(res.plaid.*repmat(dt,[1, size(res.plaid,2), size(res.plaid,3)]),1))'); 
    set(gca,'xscale','log','yscale','log');
    shading flat
end
if nargin==2
    if length(ind)==1
        if ind==1
            dn = diff(res.nums); dn(end+1) = dn(end);             
            semilogx(res.autotime,squeeze(sum(sum(res.plaid.*repmat(dn',[size(res.plaid,1), 1, size(res.plaid,3)]),2),3))); 
            xlabel('time [s]');
            ylabel('autocorrelation');
        elseif ind==2
            dt = diff(res.autotime); dt(end+1) = dt(end);
            loglog(res.nums,squeeze(sum(sum(res.plaid.*repmat(dt,[1, size(res.plaid,2), size(res.plaid,3)]),3),1))); 
            xlabel('number of photons');
            ylabel('frequency');
        elseif ind==3
            dn = diff(res.nums); dn(end+1) = dn(end);             
            dt = diff(res.autotime); dt(end+1) = dt(end);
            plot(res.bintau,squeeze(sum(sum(res.plaid.*repmat(dt*dn',[1, 1, size(res.plaid,3)]),1),2)))            
            xlabel('average lifetime [ns]');
            ylabel('frequency');
        end
    else 
        if ind(1)==2 & ind(2)==3
            dt = diff(res.autotime); dt(end+1) = dt(end);
            pcolor(res.nums,res.bintau,squeeze(sum(res.plaid.*repmat(dt,[1, size(res.plaid,2), size(res.plaid,3)]),1))'); 
            set(gca,'xscale','log')
            shading flat
            xlabel('number of photons')
            ylabel('average lifetime [ns]')             
        elseif ind(1)==1 & ind(2)==3                
            dn = diff(res.nums); dn(end+1) = dn(end);             
            %pcolor(res.autotime,res.bintau,squeeze(sum(res.plaid.*repmat(dn',[size(res.plaid,1), 1, size(res.plaid,3)]),2))'); 
            tmp = squeeze(sum(sum(res.plaid.*repmat(dn',[size(res.plaid,1), 1, size(res.plaid,3)]),2),3));
            tmp = tmp-min(tmp);
            pcolor(res.autotime,res.bintau,(tmp*ones(1,length(res.bintau)))'.*squeeze(sum(res.plaid.*repmat(res.nums'.*dn',[size(res.plaid,1), 1, size(res.plaid,3)]),2))'); 
            set(gca,'xscale','log')
            shading flat
            xlabel('lag time [s]')
            ylabel('average lifetime [ns]')             
        elseif ind(1)==1 & ind(2)==2            
            dn = diff(res.nums); dn(end+1) = dn(end);                         
            tmp = squeeze(sum(sum(res.plaid.*repmat(dn',[size(res.plaid,1), 1, size(res.plaid,3)]),2),3));
            tmp = tmp-min(tmp);
            pcolor(res.autotime,res.nums,log10((tmp*res.nums')'.*squeeze(sum(res.plaid,3))')); 
            set(gca,'xscale','log','yscale','log')
            shading flat
            xlabel('lag time [s]')
            ylabel('number of photons')             
        end
    end
end


