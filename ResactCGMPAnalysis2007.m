load c:\Joerg\Doc\Kaupp\Bungert\IngoNewdata\WeyandtData2007-10.mat
clear final
bild = false;

if 0
    for cnt = 1:length(data)
        clear t y c p1 p2 z1 z2 sy

        t = data(cnt).t';
        conc = unique(data(cnt).conc);
        clear y
        for j=1:length(conc)
            tmp = data(cnt).x(:,data(cnt).conc==conc(j));
            ind = isnan(tmp);
            tmp(ind) = 0;
            y(:,j) = sum(tmp,2)./sum(double(~ind),2);
            sy(:,j) = sqrt(sum((tmp>0).*(tmp-y(:,j)*ones(1,size(tmp,2))).^2,2)./(sum(double(~ind),2)-1));
        end
        sy(isnan(sy)) = 0;

        for j=1:size(y,2)
            ind = ~isnan(y(:,j));
            p1(j) = 0.2;
            for sub = 1:3
                p1(j) = simplex('ExpFun',p1(j),[],[],[],[],t(ind),y(ind,j));
            end
            [tmp tmp] = ExpFun(p1(:,j),t(ind),y(ind,j));
            z1(:,j) = ExpFun(p1(:,j),tmp,t);
        end

        if 0
            errorbar(t*ones(1,size(y,2)),y,sy,'o','linewidth',1.5); hold on; plot(t,z1); hold off
            ax = axis; ax(3) = 0; ax(2) = 1.8; axis(ax)
            %     for j=1:size(y,2)
            %         text(1.65,z1(end,j),['\itk_i\rm = ' mnum2str(1/p1(j),2,2) '\cdots^{-1}'])
            %     end
            s = [num2str(1e9*conc',3) repmat(' nM',length(conc),1)];
            h = legend(s,2); set(h,'fontsize',12)
            xlabel('time [s]')
            ylabel('cGMP [pmol / 10^8 sperm]');

            eval(['print -dpng -r300 res' data(cnt).name 'a'])
        end

        for j=1:size(y,2)
            ind = ~isnan(y(:,j));
            p2(:,j) = [0.1 0.2];
            for sub = 1:25
                p2(:,j) = simplex('BenFun',p2(:,j),[],[0.5 0.5],[],[],t(ind),y(ind,j));
            end
            [tmp tmp] = BenFun(p2(:,j),t(ind),y(ind,j));
            z2(:,j) = BenFun(p2(:,j),tmp,t,1); drawnow
        end
        p2 = sort(p2);

        if 1
            errorbar(t*ones(1,size(y,2)),y,sy,'o','linewidth',1.5); hold on; plot(t,z2); hold off
            ax = axis; ax(3) = 0; ax(2) = 1.8; axis(ax)
            %     for j=1:size(y,2)
            %         text(1.65,z2(end,j),{['\itk_a\rm = ' mnum2str(1/p2(1,j),2,2) '\cdots^{-1}'],...
            %             ['\itk_i\rm = ' mnum2str(1/p2(2,j),2,2) '\cdots^{-1}']})
            %     end
            s = [num2str(1e9*conc',3) repmat(' nM',length(conc),1)];
            h = legend(s,2); set(h,'fontsize',12)
            xlabel('time [s]')
            ylabel('cGMP [pmol / 10^8 sperm]');

            %eval(['print -dpng -r300 res' data(cnt).name 'b'])
        end

        final(cnt).name = data(cnt).name;
        final(cnt).conc = conc;
        final(cnt).sperm = data(cnt).sperm;
        final(cnt).tau1 = p1;
        final(cnt).tau2 = p2;
        final(cnt).p1 = 1./p1;
        final(cnt).p2 = 1./p2;
        final(cnt).c2 = tmp;

    end
    save c:\Joerg\Doc\Kaupp\Bungert\IngoNewdata\ResactCGMPResults2007 final 
end

if 0
    load c:\Joerg\Doc\Kaupp\Bungert\IngoNewdata\ResactCGMPResults2007 
    conc = []; pp = [];
    for cnt = 1:length(final)
        conc = [conc final(cnt).conc];
        pp = [pp final(cnt).p2(2,:)];
    end
    conc = round(1e12*conc);
    uc = unique(conc);
    clear ppav ppstd
    for j=1:length(uc)
        ppav(j) = mean(pp(conc==uc(j)));
        ppstd(j) = std(pp(conc==uc(j)));
    end
    c = polyfit(log(uc),log(ppav),1);
    semilogy(log10(1e-3*uc),ppav,'o',log10(1e-3*uc),exp(polyval(c,log(uc))))
    box off
    hold on
    for j=1:length(uc)
        plot(log10(1e-3*uc(j)*[1 1]), ppav(j)+[-1 1]*ppstd(j))
        plot(log10(1e-3*uc(j)*[1-0.1 1+0.1]), ppav(j)+[-1 -1]*ppstd(j))
        plot(log10(1e-3*uc(j)*[1-0.1 1+0.1]), ppav(j)+[1 1]*ppstd(j))
    end
    hold off
    xlabel('log [resact] (nM)')
    ylabel('inactivation rate (1/s)')
end

if 1
    load c:\Joerg\Doc\Kaupp\Bungert\IngoNewdata\ResactCGMPResults2007 
    NatConst
    conc = []; uu = []; pa = [];
    for cnt = 1:length(data)
        conc = [conc final(cnt).conc];
        uu = [uu final(cnt).conc/final(cnt).sperm*AvogadroConstant/1e3];
        pa = [pa final(cnt).c2(2,:).*diff(final(cnt).p2)./final(cnt).p2(1,:)*AvogadroConstant*1e-20];
    end
    uc = unique(uu);
    for j=1:length(uc)
        paav(j) = mean(pa(uu==uc(j)))/uc(j);
        pastd(j) = std(pa(uu==uc(j)))/uc(j);
    end

    c = polyfit(log(uc),log(paav),1);
    semilogy(log10(uc),paav,'o',log10(uc),exp(polyval(c,log(uc))));
    box off
    hold on
    %     for j=1:length(uc)
    %         semilogy(log10(uc(j)*[1 1]), paav(j)+[-1 1]*pastd(j))
    %         semilogy(log10(uc(j)*[1-0.1 1+0.1]), paav(j)+[-1 -1]*pastd(j))
    %         semilogy(log10(uc(j)*[1-0.1 1+0.1]), paav(j)+[1 1]*pastd(j))
    %     end
    ax = axis;
    plot(6*ones(1,20),linspace(ax(3),ax(4),20),'--k','linewidth',0.5)
    hold off
    xlabel('log [resact] (molecules/sperm)')
    ylabel('\itk_s\rm\cdotG^°/resact (1/s)')
end

