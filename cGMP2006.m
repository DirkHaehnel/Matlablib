clear all
close all
path = 'D:\Joerg\Doc\Kaupp\Bungert\cGMP IBMX Werte aus 2006\';
names = dir([path '*IBMX.mat']);

if 0
    bild = 0;
    for cnt = 1:length(names)
        clear t y c p2 
        s = int2str(cnt);

        load([path names(cnt).name])
        x(:,conc==0) = [];
        conc(conc==0) = [];
        c = conc;
        conc = union(conc,conc);
        t = time;
        clear y
        for j=1:length(conc)
            tmp = x(:,c==conc(j));
            ind = isnan(tmp);
            tmp(ind) = 0;
            y(:,j) = sum(tmp,2)./sum(double(~ind),2);
            sy(:,j) = sqrt(sum((tmp>0).*(tmp-y(:,j)*ones(1,size(tmp,2))).^2,2)./(sum(double(~ind),2)-1));
        end
        sy(isnan(sy)) = 0;
        ind = t<=t(end);
        for j=1:size(y,2)
            %ind = ~isnan(y(:,j));
            p2(:,j) = [0.5 0.2];
            for sub = 1:25
                p2(:,j) = simplex('BenFun',p2(:,j),[],[],[],[],t(ind),y(ind,j));
            end
            [tmp tmp] = BenFun(p2(:,j),t(ind),y(ind,j));
            z2(:,j) = BenFun(p2(:,j),tmp,linspace(t(1),t(end),100));
        end
        p2 = sort(p2);

        if bild
            plot(linspace(t(1),t(end),100),z2); hold on; errorbar(t*ones(1,size(y,2)),y,sy,'o','linewidth',1.5); hold off
            ax = axis; ax(1) = 0; ax(3) = 0; ax(2) = 1.62; axis(ax)
            set(gca,'tickdir','out')
            box off
            s = [num2str(conc',3) repmat(' nM',length(conc),1)];
            h = legend(s,2); set(h,'fontsize',12)
            xlabel('time (s)')
            ylabel('cGMP (pmol / 10^8 sperm)');
            %pause
            %eval(['print -dpng res' name{cnt}])
        end
    end
end

if 0
    bild = 1;
    for cnt = 1:length(names)
        clear pp
        s = int2str(cnt);

        load([path names(cnt).name])
        x(:,conc==0) = [];
        conc(conc==0) = [];
        c = conc;
        conc = union(conc,conc);
        t = time;

        for j=1:size(x,2)
            ind = ~isnan(x(:,j));
            pp(:,j) = [0.1 0.2];
            for sub = 1:15
                pp(:,j) = simplex('BenFun',pp(:,j),[],[],[],[],t(ind),x(ind,j));
            end
            pp(:,j) = sort(pp(:,j));
            [tmp cc(:,j)] = BenFun(pp(:,j),t(ind),x(ind,j),bild);
        end

        final(cnt).pp = 1./pp; % BenFun retruns dwell time, not rate !!!
        final(cnt).cc = cc;
        final(cnt).conc = c;
        final(cnt).t = t;
    end
    eval(['save ''' path '''Results final']);
end
    
if 1
    load([path 'Results'])
    conc = []; pp = [];
    for cnt = 1:4%length(final)
        conc = [conc final(cnt).conc];
        pp = [pp final(cnt).pp(2,:)];
    end
    c = polyfit(log(conc),log(pp),1);
    uc = union(conc,conc);
    clear ppav ppstd
    for j=1:length(uc)
        ppav(j) = mean(pp(conc==uc(j)));
        ppstd(j) = std(pp(conc==uc(j)));
    end
    semilogy(log10(uc),ppav,'o',log10(uc),exp(polyval(c,log(uc))))
    box off
    hold on
    for j=1:length(uc)
        plot(log10(uc(j)*[1 1]), ppav(j)+[-1 1]*ppstd(j))
        plot(log10(uc(j)*[1-0.1 1+0.1]), ppav(j)+[-1 -1]*ppstd(j))
        plot(log10(uc(j)*[1-0.1 1+0.1]), ppav(j)+[1 1]*ppstd(j))
    end
    hold off
    xlabel('log [resact] (nM)')
    ylabel('inactivation rate (1/s)')
end

if 0
    NatConst
    spermv = [0.97 1.68 0.965 0.885]*1e8;
    load([path 'Results'])
    conc = []; pp = []; uu = []; vv = [];
    for cnt = 1:4%length(final)
        conc = [conc final(cnt).conc];
        pp = [pp final(cnt).pp(2,:)];
        uu = [uu final(cnt).conc*1e-9/spermv(cnt)*AvogadroConstant/1e3]; % number of resact per sperm
        vv = [vv final(cnt).cc(2,:).*diff(final(cnt).pp)./final(cnt).pp(1,:)*AvogadroConstant*1e-20]; % y(t) was pmol cGMP/10^8 sperm
    end
    uc = union(uu,uu);
    clear paav ppstd
    for j=1:length(uc)
        paav(j) = mean(vv(uu==uc(j)))/uc(j);
        pastd(j) = std(vv(uu==uc(j)))/uc(j);
    end

    act = polyfit(log(uc),log(paav),1);
    semilogy(log10(uc),paav,'o',log10(uc),exp(polyval(act,log(uc))));
    box off
    hold on
    for j=1:length(uc)
        semilogy(log10(uc(j)*[1 1]), paav(j)+[-1 1]*pastd(j))
        semilogy(log10(uc(j)*[1-0.1 1+0.1]), paav(j)+[-1 -1]*pastd(j))
        semilogy(log10(uc(j)*[1-0.1 1+0.1]), paav(j)+[1 1]*pastd(j))
    end
    ax = axis;
    plot(6*ones(1,20),linspace(ax(3),ax(4),20),'--k','linewidth',0.5)
    hold off
    xlabel('log [resact] (molecules/sperm)')
    ylabel('\itk_s\rm\cdotG^*/resact (1/s)')
end
