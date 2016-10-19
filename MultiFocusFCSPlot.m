function MultiFocusFCSPlot(y, t)

if ((size(y,2) == 6) || (size(y,2) == 4))
    mode = 3;
    t1 = 'ACF + CCF';
elseif (size(y,2) == 16)
    mode = 2;
    t1 = 'blue ACF + CCF';
    t2 = 'FRET ACF + CCF';
    t3 = '2 Color CCF + 2 Color - 2 Focus CCCF';    
elseif ((size(y,2) == 22) || (size(y,2) == 26))
    mode = 1;
    t1 = 'red ACF + CCF';
    t2 = 'blue ACF + CCF';
    t3 = 'FRET ACF + CCF';
    t4 = '2 Color CCF + 2 Color - 2 Focus CCCF';
else
    fprintf(1,'\n\n      Wrong input size\n');
    return
end

clf;

if (mode == 3)
    if (size(y,2) == 6)
        y(:,3,:)=mean(y(:,[3 5],:),2);
        y(:,4,:)=mean(y(:,[4 6],:),2);
        y(:,5:6,:)=[];
    end
    semilogx(1e6*t,mean(y(:,1:4,:),3));
    ax=axis;
    axis([1e6*min(t) 1e6*max(t) min([min(min(mean(y(:,1:4,:),3))) min(abs(ax(3:4)))]) max(abs(ax(3:4)))]);
    ylabel('correlation / (\ith\rm\nu)^2/s^2');
    xlabel('time [\mus]');
    title(t1);
else
    %text(0.6,0.95,s,'VerticalAlignment','Top','units','normal')
    subplot('Position',[0.075 0.575 0.425 0.35]);
    y(:,3,:)=mean(y(:,[3 5],:),2);
    y(:,4,:)=mean(y(:,[4 6],:),2);
    semilogx(1e6*t,mean(y(:,1:4,:),3));
    ax=axis;
    axis([1e6*min(t) 1e6*max(t) min([min(min(mean(y(:,1:4,:),3))) min(abs(ax(3:4)))]) max(abs(ax(3:4)))]);
    ylabel('correlation / (\ith\rm\nu)^2/s^2');
    title(t1);

    subplot('Position',[0.55 0.575 0.425 0.35]);
    y(:,9,:)=mean(y(:,[9 11],:),2);
    y(:,10,:)=mean(y(:,[10 12],:),2);
    semilogx(1e6*t,mean(y(:,7:10,:),3));
    ax=axis;
    axis([1e6*min(t) 1e6*max(t) min([min(min(mean(y(:,7:10,:),3))) min(abs(ax(3:4)))]) max(abs(ax(3:4)))]);
    title(t2);

    subplot('Position',[0.075 0.1 0.425 0.35]);
    if (mode == 1)
        y(:,15,:)=mean(y(:,[15 17],:),2);
        y(:,16,:)=mean(y(:,[16 18],:),2);
    end
    semilogx(1e6*t,mean(y(:,13:16,:),3));
    ax=axis;
    axis([1e6*min(t) 1e6*max(t) min([min(min(mean(y(:,13:16,:),3))) min(abs(ax(3:4)))]) max(abs(ax(3:4)))]);
    ylabel('correlation / (\ith\rm\nu)^2/s^2');
    xlabel('time [\mus]');
    title(t3);

    if (mode == 1)
        subplot('Position',[0.55 0.1 0.425 0.35]);
        semilogx(1e6*t,mean(y(:,19:22,:),3));
        ax=axis;
        axis([1e6*min(t) 1e6*max(t) min([min(min(mean(y(:,19:22,:),3))) min(abs(ax(3:4)))]) max(abs(ax(3:4)))]);
        xlabel('time [\mus]');
        title(t4);
    end

    set(gcf,'Position',[20 50 1200 850],'PaperPositionMode','auto')
end

