load D:\Joerg\Doc\IngoW\cGMP_IBMX.mat
% load D:\Doc\IngoW\cGMP_IBMX.mat
clear final
bild = true;

if 0 
    for cnt = 1:10
        clear t y c p1 p2 z1 z2 sy
        s = int2str(cnt);

        eval(['t = t' s ';']);
        eval(['conc = union(c' s ',c' s ');']);
        clear y
        for j=1:length(conc)
            eval(['tmp = x' s '(:,c' s '==conc(j));']);
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

        if bild
            errorbar(t*ones(1,size(y,2)),y,sy,'o','linewidth',1.5); hold on; plot(t,z1); hold off
            ax = axis; ax(3) = 0; ax(2) = 1.8; axis(ax)
            %     for j=1:size(y,2)
            %         text(1.65,z1(end,j),['\itk_i\rm = ' mnum2str(1/p1(j),2,2) '\cdots^{-1}'])
            %     end
            s = [num2str(1e9*conc',3) repmat(' nM',length(conc),1)];
            h = legend(s,2); set(h,'fontsize',12)
            xlabel('time [s]')
            ylabel('cGMP [pmol / 10^8 sperm]');

            eval(['print -dpng -r300 res' name{cnt} 'a'])
        end

        for j=1:size(y,2)
            ind = ~isnan(y(:,j));
            p2(:,j) = [0.1 0.2];
            for sub = 1:25
                p2(:,j) = simplex('BenFun',p2(:,j),[],[0.5 0.5],[],[],t(ind),y(ind,j));
            end
            [tmp tmp] = BenFun(p2(:,j),t(ind),y(ind,j));
            z2(:,j) = BenFun(p2(:,j),tmp,t);
        end
        p2 = sort(p2,2);

        if bild
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

            eval(['print -dpng -r300 res' name{cnt} 'b'])
        end

        final(cnt).tau1 = p1;
        final(cnt).tau2 = p2;
        final(cnt).p1 = 1./p1;
        final(cnt).p2 = 1./p2;

    end
else
    for cnt = 1:10
        clear t y c p1 p2 z1 z2 sy
        s = int2str(cnt);

        eval(['t = t' s ';']);
        eval(['conc = c' s ';']);
        eval(['y = x' s ';']);

        for j=1:size(y,2)
            ind = ~isnan(y(:,j));
            p2(:,j) = [0.1 0.2];
            for sub = 1:15
                p2(:,j) = simplex('BenFun',p2(:,j),[],[0.5 0.5],[],[],t(ind),y(ind,j));
            end
            [tmp tmp] = BenFun(p2(:,j),t(ind),y(ind,j),bild);
            z2(:,j) = BenFun(p2(:,j),tmp,t);
        end
        p2 = sort(p2);

        final(cnt).tau2 = p2;
        final(cnt).p2 = 1./p2;
    end
    save D:\Joerg\Doc\IngoW\ResactCGMPResults_all.mat final 
end
    

return

fout=fopen('tst.dat','w'); 
fprintf(fout,'%s\n\n\n','p1');
for j=1:10 
    fprintf(fout,'%s\n\n',name{j});
    for k=1:size(final(j).p1,2)
        fprintf(fout,'%d ', final(j).p1(k));
    end;
    fprintf(fout,'\n\n');
end
fprintf(fout,'\n\n%s\n\n\n','tau1');
for j=1:10 
    fprintf(fout,'%s\n\n',name{j});
    for k=1:size(final(j).p1,2)
        fprintf(fout,'%d ', final(j).tau1(k));
    end
    fprintf(fout,'\n\n'); 
end

fprintf(fout,'\n\n%s\n\n\n','p2');
for j=1:10 
    fprintf(fout,'%s\n\n',name{j});
    for kk=1:2 
        for k=1:size(final(j).p2,2) 
            fprintf(fout,'%d ', final(j).p2(kk,k)); 
        end; 
        fprintf(fout,'\n'); 
    end
    fprintf(fout,'\n'); 
end
fprintf(fout,'\n\n%s\n\n\n','tau2');
for j=1:10 
    fprintf(fout,'%s\n\n',name{j});
    for kk=1:2 
        for k=1:size(final(j).p2,2) 
            fprintf(fout,'%d ', final(j).tau2(kk,k)); 
        end; 
        fprintf(fout,'\n'); 
    end
    fprintf(fout,'\n'); 
end
fclose(fout)

return

for cnt = 1:10
    s = int2str(cnt);

    eval(['t = t' s ';']);
    eval(['conc = union(c' s ',c' s ');']);

    if cnt==1
        loglog(conc,final(cnt).p2(2,:),'-')
        hold on
    else
        loglog(conc,final(cnt).p2(2,:),'-')
    end
end
hold off
colorize(get(gca,'children'))
xlabel('resact concentration [M]')
ylabel('inactivation rate [1/s]')
legend(name,-1)

for cnt = 1:9
    if ~(cnt==3)
        s = int2str(cnt);

        eval(['t = t' s ';']);
        eval(['conc = union(c' s ',c' s ');']);

        if cnt==1
            loglog(conc,final(cnt).p2(1,:)/spermv(cnt),'-')
            hold on
        else
            loglog(conc,final(cnt).p2(1,:)/spermv(cnt),'-')
        end
    end
end
hold off
colorize(get(gca,'children'))
xlabel('resact concentration [M]')
ylabel('activation rate [ml/s/sperms]')
legend(name{[1:2 4:9]},-1)

for cnt = 1:10
    s = int2str(cnt);

    eval(['t = t' s ';']);
    eval(['conc = union(c' s ',c' s ');']);

    if cnt==1
        loglog(conc,final(cnt).p1,'-')
        hold on
    else
        loglog(conc,final(cnt).p1,'-')
    end
end
hold off
colorize(get(gca,'children'))
xlabel('resact concentration [M]')
ylabel('inactivation rate [1/s]')
legend(name,-1)


return % Kaupp Images

load D:\Joerg\Doc\IngoW\cGMP_IBMX.mat
load D:\Joerg\Doc\IngoW\ResactCGMPResults_all.mat

conc = []; p2 = [];
for cnt = 3:10 % omit Johannes' data 2000
    s = int2str(cnt);
    %eval(['conc = [conc union(c' s ',c' s ')];']);
    eval(['conc = [conc c' s '];']);
	p2 = [p2 final(cnt).p2(2,:)];
end
uc = union(conc,conc);
for j=1:length(uc)
    p2av(j) = mean(p2(conc==uc(j)));
    p2std(j) = std(p2(conc==uc(j)));
end
c = polyfit(log(uc),log(p2av),1);
loglog(uc,p2av,'o',uc,exp(polyval(c,log(uc))))
hold on
for j=1:length(uc)
    plot(uc(j)*[1 1], p2av(j)+[-1 1]*p2std(j))
    plot(uc(j)*[1-0.1 1+0.1], p2av(j)+[-1 -1]*p2std(j))
    plot(uc(j)*[1-0.1 1+0.1], p2av(j)+[1 1]*p2std(j))    
end
hold off
xlabel('resact concentration [M]')
ylabel('inactivation rate [1/s]')



for cnt = 4:9
    s = int2str(cnt);
    eval(['t = t' s ';']);
    eval(['conc = c' s ';']);
    p2 = [p2 final(cnt).p2(1,:)/spermv(cnt)];
end
hold off
colorize(get(gca,'children'))
xlabel('resact concentration [M]')
ylabel('activation rate [ml/s/sperms]')
legend(name{[1:2 4:9]},-1)
