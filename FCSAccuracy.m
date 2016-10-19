% Program for calculating the accuracy of the autocorrelation function

if 0
    
    a = 0.5; % small half axis of MDF ellipsoid
    b = 2; % large half axis of MDF ellipsoid
    N = 200; % max number of time channels
    dt = 1e-8; % autocorrelation time interval in seconds
    tmax = 1; % max autocorrelation time in seconds
    kappa = log(tmax/dt)/(N-1);
    t = dt*exp((0:N-1)*kappa);
    Nc = 100; cmin = 1e-5; cmax = 1e3;
    kappa = log(cmax/cmin)/(Nc-1);
    cv = cmin*exp((0:Nc-1)*kappa); % concentration vector in molecules per mum^3
    kappav = 10.^(1:0.1:5); % photons per second
    meastime = 1e2:1e2:1e4; % measurement time in seconds
    
    % diffusion constants to be considered
    dmin = 1e2;
    dif0 = 1e3;
    nd = 100;
    kappa = log(dif0/dmin)/nd*2;
    dif = dmin*exp(kappa*(0:nd)); dif = dif'; % diffusion constant in mum^2/second unit
    
    Z0 = FCSZfun(a,b,dif0);
    
    for jc = 1:length(cv)
        dif(nd/2+1) = [];       
        c = cv(jc); 
        
        [w, w0] = FCSalpha(dif0,dif,c,a,b,t);
        %    alpha = FCSalpha(dif0,dif,c,a,b,t,w0*t',w*t');
        alpha = FCSalpha(dif0,dif,c,a,b,t,sum(w0),sum(w')')./(ones(size(dif))*t);
        
        M = (c*FCSZfun(a,b,dif0,t) + c^2*Z0^2).*t;
        MM2 = alpha.^2*M';
        M = alpha*M';
        
        MM3 = 0*dif(:);
        for u=1:N
            for v=1:N
                t1 = t(u); t2 = t(v);        
                fac = t1*t2*alpha(:,u).*alpha(:,v);        
                MM3 = MM3 + fac*(c*FCSZfun(a,b,dif0,t1,t2)+...
                    c^2*Z0*(FCSZfun(a,b,dif0,t1)+FCSZfun(a,b,dif0,t2)+FCSZfun(a,b,dif0,t1+t2))+...
                    c^3*Z0^3);
                if u<v
                    t2 = t(v) - t1;
                    MM3 = MM3 + 2*fac*(c*FCSZfun(a,b,dif0,t1,t2)+...
                        c^2*Z0*(FCSZfun(a,b,dif0,t1)+FCSZfun(a,b,dif0,t2)+FCSZfun(a,b,dif0,t1+t2))+...
                        c^3*Z0^3);
                end
            end
            u
        end
        MM3 = 2*MM3;
        
        MM4 = 0*dif(:);
        for u=1:N
            for w=1:N
                t1 = t(u); t3 = t(w); 
                fac = t1*t3*alpha(:,u).*alpha(:,w); 
                t2 = t; tt = t';
                MM4 = MM4 + fac*(c*FCSZfun(a,b,dif0,t1,t2,t3)*tt+...
                    c^2*Z0*(FCSZfun(a,b,dif0,t1,t2)+FCSZfun(a,b,dif0,t1,t2+t3)+...
                    FCSZfun(a,b,dif0,t1+t2,t3)+FCSZfun(a,b,dif0,t2,t3))*tt+...
                    c^2*(FCSZfun(a,b,dif0,t1+t2).*FCSZfun(a,b,dif0,t2+t3)+...
                    FCSZfun(a,b,dif0,t1+t2+t3).*FCSZfun(a,b,dif0,t2))*tt+...
                    c^3*Z0^2*(FCSZfun(a,b,dif0,t2)+...
                    FCSZfun(a,b,dif0,t1+t2)+FCSZfun(a,b,dif0,t2+t3)+FCSZfun(a,b,dif0,t1+t2+t3))*tt);
                ind = 1:N;
                ind = ind(u>ind & w>ind);            
                if ~isempty(ind)
                    tt = t(ind)';
                    t1 = t(u) - t(ind);
                    t2 = t(ind);
                    t3 = t(w) - t(ind);
                    MM4 = MM4 + fac*(c*FCSZfun(a,b,dif0,t1,t2,t3)*tt+...
                        c^2*Z0*(FCSZfun(a,b,dif0,t1,t2)+FCSZfun(a,b,dif0,t1,t2+t3)+...
                        FCSZfun(a,b,dif0,t1+t2,t3)+FCSZfun(a,b,dif0,t2,t3))*tt+...
                        c^2*(FCSZfun(a,b,dif0,t1).*FCSZfun(a,b,dif0,t3)+FCSZfun(a,b,dif0,t1+t2).*FCSZfun(a,b,dif0,t2+t3)+...
                        FCSZfun(a,b,dif0,t1+t2+t3).*FCSZfun(a,b,dif0,t2))*tt+...
                        c^3*Z0^2*(FCSZfun(a,b,dif0,t1)+FCSZfun(a,b,dif0,t2)+FCSZfun(a,b,dif0,t3)+...
                        FCSZfun(a,b,dif0,t1+t2)+FCSZfun(a,b,dif0,t2+t3)+FCSZfun(a,b,dif0,t1+t2+t3))*tt+...
                        c^4*Z0^4*sum(tt));
                end
                ind = 1:N;
                ind = ind(u-w>ind);
                if ~isempty(ind)
                    tt = t(ind)';
                    t1 = t(u) - t(ind) - t(w);
                    t2 = t(w);
                    t3 = t(ind);
                    MM4 = MM4 + fac*(c*FCSZfun(a,b,dif0,t1,t2,t3)*tt+...
                        c^2*Z0*(FCSZfun(a,b,dif0,t1,t2)+FCSZfun(a,b,dif0,t1,t2+t3)+...
                        FCSZfun(a,b,dif0,t1+t2,t3)+FCSZfun(a,b,dif0,t2,t3))*tt+...
                        c^2*(FCSZfun(a,b,dif0,t1).*FCSZfun(a,b,dif0,t3)+FCSZfun(a,b,dif0,t1+t2).*FCSZfun(a,b,dif0,t2+t3)+...
                        FCSZfun(a,b,dif0,t1+t2+t3).*FCSZfun(a,b,dif0,t2))*tt+...
                        c^3*Z0^2.*(FCSZfun(a,b,dif0,t1)+FCSZfun(a,b,dif0,t2)+FCSZfun(a,b,dif0,t3)+...
                        FCSZfun(a,b,dif0,t1+t2)+FCSZfun(a,b,dif0,t2+t3)+FCSZfun(a,b,dif0,t1+t2+t3))*tt+...
                        c^4*Z0^4*sum(tt));
                end        
            end
            u
        end        
        
        err=0*dif;
        for jk = 1:length(kappav)
            kappa = kappav(jk);
            for jt = 1:length(meastime)
                tt = meastime(jt);
                dM = (kappa^2*MM2 + kappa^3*MM3 + kappa^4*MM4)/tt;
                err(:,jk,jt) = 1 - erf(kappa^2*M./sqrt(2*dM));
            end
        end
        err = cat(1,err,ones(size(err(1,:,:))));
        err = err([1:nd/2 end nd/2+1:nd],:,:);
        dif = [dif(1:nd/2); dif0; dif(nd/2+1:100)];
        
        save(['fcsaccures' mint2str(jc,3)],'M','MM2','MM3','MM4','err','kappav','meastime','c','dif0','nd','dif');
    end
    
end

if 0
    % M-value figure
    a = 0.5; % small half axis of MDF ellipsoid
    b = 2; % large half axis of MDF ellipsoid
    N = 200; % max number of time channels
    dt = 1e-8; % autocorrelation time interval in seconds
    tmax = 1; % max autocorrelation time in seconds
    kappa = log(tmax/dt)/(N-1);
    t = dt*exp((0:N-1)*kappa);
    Nc = 100; cmin = 1e-5; cmax = 1e3;
    kappa = log(cmax/cmin)/(Nc-1);
    cv = cmin*exp((0:Nc-1)*kappa); % concentration vector in molecules per mum^3
    kappav = 10.^(1:0.1:5); % photons per second
    meastime = 1e2:1e2:1e4; % measurement time in seconds
    
    % diffusion constants to be considered
    dmin = 1e2;
    dif0 = 1e3;
    nd = 100;
    kappa = log(dif0/dmin)/nd*2;
    dif = dmin*exp(kappa*(0:nd)); dif = dif'; % diffusion constant in mum^2/second unit
    
    Z0 = FCSZfun(a,b,dif0);

    jc = 63;
    dif(nd/2+1) = [];       
    c = cv(jc); 
    
    [w, w0] = FCSalpha(dif0,dif,c,a,b,t);
    %    alpha = FCSalpha(dif0,dif,c,a,b,t,w0*t',w*t');
    alpha = FCSalpha(dif0,dif,c,a,b,t,sum(w0),sum(w')')./(ones(size(dif))*t);
    
    M = (c*FCSZfun(a,b,dif0,t) + c^2*Z0^2).*t;
    MM2 = alpha.^2*M';
    M = alpha*M';
    
    semilogx(dif,(alpha(10:10:40,:).*(ones(4,1)*t))*w',dif,0*dif,'k:')
    hold on; 
    plot([dif(10) dif(10)],[-0.5 0.5],':r'); 
    plot([dif(20) dif(20)],[-0.5 0.5],':b'); 
    plot([dif(30) dif(30)],[-0.5 0.5],':g'); 
    plot([dif(40) dif(40)],[-0.5 0.5],':c'); 
    plot([dif0 dif0],[-0.5 0.5],':k'); 
    hold off
    
    xlabel('diffusion coefficient [\mum^2/s]')
    ylabel('average \itM\rm-value')

    % error distribution figure
    load fcsaccures063
    len = 8;
    semilogx(dif,squeeze(err(:,21,5)),'color',[1-0/len 0 0/len]);
    for j=2:9
        line(dif,squeeze(err(:,21,j*10-5)),'color',[1-(j-1)/len 0 (j-1)/len],'linewidth',2);
    end
    line(dif,0.1*(dif>0),'color','k','linestyle',':')
    legend(num2str(meastime(5:10:end-10)'))
    xlabel('diffusion coefficient [\mum^2/s]')
    ylabel('error')

    % width of error distribution
    for j=1:size(err,3)
        for k=1:size(err,2)
            tmp=dif(err(:,k,j)>0.1); 
            res1(j,k) = (tmp(end)/dif0+dif0/tmp(1))/2;
        end
    end
    clear c1
    for j=1:41; 
        c1(j,:)=polyfit([meastime(res1(:,j)<9) 2e3+2*max(meastime(res1(:,j)<9))-fliplr(meastime(res1(:,j)<9))],log(log([res1(res1(:,j)<9,j);flipud(res1(res1(:,j)<9,j))]))',8); 
    end
    len = 3;
    j=15; plot(meastime(res1(:,j)<9),exp(exp(polyval(c1(j,:),meastime(res1(:,j)<9)))),'color',[1-0/len 0 0/len])
    j=18; line(meastime(res1(:,j)<9),exp(exp(polyval(c1(j,:),meastime(res1(:,j)<9)))),'color',[1-1/len 0 1/len])
    j=21; line(meastime(res1(:,j)<9),exp(exp(polyval(c1(j,:),meastime(res1(:,j)<9)))),'color',[1-2/len 0 2/len])
    j=41; line(meastime(res1(:,j)<9),exp(exp(polyval(c1(j,:),meastime(res1(:,j)<9)))),'color',[1-3/len 0 3/len])
    axis([0 1e4 1 5])    
    xlabel('measurement time [s]'); ylabel('0.5\times(\itD_{\rmmax}/\itD_{\rm0}+\itD_{\rm0}/\itD_{\rmmin}\rm)');
    legend({'250','500','10^3','10^4'})
end

if 1
    a = 0.5; % small half axis of MDF ellipsoid
    b = 2; % large half axis of MDF ellipsoid
    V = pi^(3/2)*a^2*b;    
    load fcsaccures001
    k1 = sum(dif<dif0/2)+1;
    k2 = sum(dif<2*dif0);
    clear cerr1 cerr2
    for j=1:100
        eval(['load fcsaccures' mint2str(j,3)]);
        conc(j) = c;
        cerr1(j,:,:) = err(k1,:,:);
        cerr2(j,:,:) = err(k2,:,:);
        j
    end
    len = length(kappav);
    c1 = {}; c2 = {}; tsh = 0.1; %conc = conc*V;
    for j=1:len
        c1{j} = contour(conc,meastime,squeeze(cerr1(:,j,:))',tsh*[1 1]);
        c2{j} = contour(conc,meastime,squeeze(cerr2(:,j,:))',tsh*[1 1]);
    end
    
    ind = [18 21 24 28 31 34 38 41]
    len = length(ind)+1;
    cnt = 1; j=1;
    semilogx(c1{ind(j)}(1,cnt+(1:c1{ind(j)}(2,cnt))),c1{ind(j)}(2,cnt+(1:c1{ind(j)}(2,cnt))),'color',[1-(j-1)/len 0 (j-1)/len],'linewidth',2);
    cnt = cnt + c1{ind(j)}(2,cnt) +1;
    while cnt<=size(c1{ind(j)},2)
        line(c1{ind(j)}(1,cnt+(1:c1{ind(j)}(2,cnt))),c1{ind(j)}(2,cnt+(1:c1{ind(j)}(2,cnt))),'color',[1-(j-1)/len 0 (j-1)/len],'linewidth',2);
        cnt = cnt + c1{ind(j)}(2,cnt) +1;
    end
    for j=2:length(ind)
        cnt = 1; 
        while cnt<=size(c1{ind(j)},2)
            line(c1{ind(j)}(1,cnt+(1:c1{ind(j)}(2,cnt))),c1{ind(j)}(2,cnt+(1:c1{ind(j)}(2,cnt))),'color',[1-(j-1)/len 0 (j-1)/len],'linewidth',2);
            cnt = cnt + c1{ind(j)}(2,cnt) +1;
        end
    end
    line([1e-5 1e1],[100 100],[10 10],'color',[0 0 0],'linewidth',2,'linestyle',':');
    ax = axis;
    ax(1) = 1e-5;
    ax(2) = 1e1;
    axis(ax)
    clear tmp
    % for j=1:length(ind) tmp{j}=mint2str(round(kappav(ind(j))/10^(floor(log10(kappav(ind(j))))))*10^(floor(log10(kappav(ind(j)))))); end    
    for j=1:length(ind) tmp{j}=mint2str(kappav(ind(j))); end        
    legend({'5\times10^2','10^3','2\times10^3','5\times10^3','10^4','2\times10^4','5\times10^4','10^5'},-1)
    xlabel('concentration [1/\mum^3]'); ylabel('measurement time [s]');

    
    ind = [18 21 24 28 31 34 38 41]
    len = length(ind);
    cnt = 1; j=1;
    while cnt<=size(c2{ind(j)},2)
        semilogx(c2{ind(j)}(1,cnt+(1:c2{ind(j)}(2,cnt))),c2{ind(j)}(2,cnt+(1:c2{ind(j)}(2,cnt))),'color',[1-j/len 0 j/len],'linewidth',2);
        cnt = cnt + c2{ind(j)}(2,cnt) +1;
    end
    for j=2:length(ind)
        cnt = 1; 
        while cnt<=size(c2{ind(j)},2)
            line(c2{ind(j)}(1,cnt+(1:c2{ind(j)}(2,cnt))),c2{ind(j)}(2,cnt+(1:c2{ind(j)}(2,cnt))),'color',[1-j/len 0 j/len],'linewidth',2);
            cnt = cnt + c2{ind(j)}(2,cnt) +1;
        end
    end
    line([3e-5 1e1],[100 100],[10 10],'color',[0 0 0],'linewidth',2,'linestyle',':');
    ax = axis;
    ax(1) = 3e-5;
    axis(ax)
    clear tmp
    for j=1:length(ind) tmp{j}=mint2str(round(kappav(ind(j))/10^(floor(log10(kappav(ind(j))))))*10^(floor(log10(kappav(ind(j)))))); end    
    legend({'5\times10^2','10^3','2\times10^3','5\times10^3','10^4','2\times10^4','5\times10^4','10^5'},-1);
    xlabel('concentration [1/detection volume]'); ylabel('measurement time [s]');

end