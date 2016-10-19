if 1 % simple OTF
    
    [q,w] = meshgrid(0.0025:0.005:4,-4:0.005:4);
    amp0 = zeros(size(q));
    amp0(q.^2+w.^2>0.98 & q.^2+w.^2<1.02 & w<0 & abs(atan(q./w))<asin(1.2/1.33))=1;
    
    qq = 0:0.5:2e3; nmax = numel(qq);
    qv = q(1,:)';
    btrans = zeros(numel(qv),nmax);
    for j=1:nmax
        btrans(:,j) = besselj(0,qv/max(qv)*qq(j));
    end
    
    ff = amp0*(btrans.*repmat(qv,[1 nmax]));
    ff = mConv(ff,flipud(ff));
    amp1 = ff*(btrans.*repmat(qq,[numel(qv) 1]))';
    
    amp1 = amp1./max(amp1(:));
    amp1(amp1<1e-6) = 0;
    
%     CombineImages(log10(cat(3,[fliplr(amp0) amp0],[fliplr(amp1) amp1])),1,2,[],{{'(A)',''} ,{'(B)',''}})
%     h = get(gca,'children');
%     for j=1:2
%         set(h(j),'fontname','calibri','fontsize',12)
%     end
%     hold
%     plot(size(q,2), 0.5*size(q,1),'q','Markersize',10,'color','yellow','linewidth',1)
%     hold
%     line(size(q,2)*[1 1], [0.9 0.8]*size(q,1),'color','yellow','linewidth',1)
%     line(2*[0.5 0.6]*size(q,2),0.9*size(q,1)*[1 1],'color','yellow','linewidth',1)
%     text(2*0.63*size(q,2),0.9*size(q,1),'k_x','color','yellow','fontname','calibri','fontsize',10,'fontangle','italic')
%     text(2*0.5*size(q,2),0.73*size(q,1),'k_z','color','yellow','fontname','calibri','fontsize',10,'fontangle','italic','HorizontalAlignment','center')
%     caxis([-4 0])
%     h = colorbar;
%     set(h,'LineWidth',0.5)
%     set(h,'Fontname','calibri','FontSize',10)

    CombineImages(log10(cat(3,[fliplr(amp0) amp0],[fliplr(amp0) amp0],[fliplr(amp1) amp1])),1,3,[],{{'(A)',''} ,{'(B)',''},{'(C)',''}})
    h = get(gca,'children');
    for j=1:3
        set(h(j),'fontname','calibri','fontsize',8)
    end
    hold
    plot(size(q,2), 0.5*size(q,1),'x','Markersize',7,'color','yellow','linewidth',1)
    hold
    line(2*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    line(4*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    hold on
    plot(2*size(q,2),size(q,1)*0.5,'xy','markersize',15,'linewidth',1); 
    plot(2*size(q,2),size(q,1)*0.5,'oy','markersize',15,'linewidth',1); 
    hold off
    text(4*size(q,2),0.47*size(q,1),'=','color','yellow','fontname','symbol','fontsize',24,'horizontalalignment','center')
    line(size(q,2)*[1 1]*0.2, [0.9 0.8]*size(q,1),'color','yellow','linewidth',1)
    line([0.2 0.4]*size(q,2),0.9*size(q,1)*[1 1],'color','yellow','linewidth',1)
    text(0.46*size(q,2),0.9*size(q,1),'k_x','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic')
    text(0.2*size(q,2),0.73*size(q,1),'k_z','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic','HorizontalAlignment','center')
    caxis([-4 0])
    h = colorbar;
    set(h,'LineWidth',0.5)
    set(h,'Fontname','calibri','FontSize',7)
    set(h,'PlotBoxAspectRatio',[1 10 1])
    
    print -r600 -dpng OTFWidefield

end

if 1 % confocal OTF
    
    [q,w] = meshgrid(0.0025:0.005:4,-4:0.005:4);
    amp0 = zeros(size(q));
    amp0(q.^2+w.^2>0.98 & q.^2+w.^2<1.02 & w<0 & abs(atan(q./w))<asin(1.2/1.33))=1;
    
    qq = 0:0.5:2e3; nmax = numel(qq);
    qv = q(1,:)';
    btrans = zeros(numel(qv),nmax);
    for j=1:nmax
        btrans(:,j) = besselj(0,qv/max(qv)*qq(j));
    end
    
    ff = amp0*(btrans.*repmat(qv,[1 nmax]));
    ff = mConv(ff,flipud(ff));
    amp1 = ff*(btrans.*repmat(qq,[numel(qv) 1]))';
    
    amp1 = amp1./max(amp1(:));
    amp1(amp1<1e-6) = 0;
    
    ff = amp1*(btrans.*repmat(qv,[1 nmax]));
    ff = mConv(ff,ff);
    amp2 = ff*(btrans.*repmat(qq,[numel(qv) 1]))';
    
    amp2 = amp2./max(amp2(:));
    amp2(amp2<1e-6) = 0;
    
    CombineImages(log10(cat(3,[fliplr(amp1) amp1],[fliplr(amp1) amp1],[fliplr(amp2) amp2])),1,3,[],{{'(A)',''} ,{'(B)',''},{'(C)',''}})
    h = get(gca,'children');
    for j=1:3
        set(h(j),'fontname','calibri','fontsize',8)
    end
    line(2*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    line(4*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    hold on
    plot(2*size(q,2),size(q,1)*0.5,'xy','markersize',15,'linewidth',1); 
    plot(2*size(q,2),size(q,1)*0.5,'oy','markersize',15,'linewidth',1); 
    hold off
    text(4*size(q,2),0.47*size(q,1),'=','color','yellow','fontname','symbol','fontsize',24,'horizontalalignment','center')
    line(size(q,2)*[1 1]*0.2, [0.9 0.8]*size(q,1),'color','yellow','linewidth',1)
    line([0.2 0.4]*size(q,2),0.9*size(q,1)*[1 1],'color','yellow','linewidth',1)
    text(0.46*size(q,2),0.9*size(q,1),'k_x','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic')
    text(0.2*size(q,2),0.73*size(q,1),'k_z','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic','HorizontalAlignment','center')
    caxis([-4 0])
    h = colorbar;
    set(h,'LineWidth',0.5)
    set(h,'Fontname','calibri','FontSize',7)
    set(h,'PlotBoxAspectRatio',[1 10 1])
    
    print -r600 -dpng OTFConfocal

    
end

if 1 
    
    [q,w] = meshgrid(0.0025:0.005:4,-4:0.005:4);
    amp0 = zeros(size(q));
    amp0(q.^2+w.^2>0.98 & q.^2+w.^2<1.02 & abs(atan(q./w))<asin(1.2/1.33))=1;
    
    qq = 0:0.5:2e3; nmax = numel(qq);
    qv = q(1,:)';
    btrans = zeros(numel(qv),nmax);
    for j=1:nmax
        btrans(:,j) = besselj(0,qv/max(qv)*qq(j));
    end
    
    ff = amp0*(btrans.*repmat(qv,[1 nmax]));
    ff = mConv(ff,flipud(ff));
    amp1 = ff*(btrans.*repmat(qq,[numel(qv) 1]))';
    
    amp1 = amp1./max(amp1(:));
    amp1(amp1<1e-6) = 0;
    
    CombineImages(log10(cat(3,[fliplr(amp0) amp0],[fliplr(amp0) amp0],[fliplr(amp1) amp1])),1,3,[],{{'(A)',''} ,{'(B)',''},{'(C)',''}})
    h = get(gca,'children');
    for j=1:3
        set(h(j),'fontname','calibri','fontsize',8)
    end
    line(2*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    line(4*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    hold on
    plot(2*size(q,2),size(q,1)*0.5,'xy','markersize',15,'linewidth',1); 
    plot(2*size(q,2),size(q,1)*0.5,'oy','markersize',15,'linewidth',1); 
    hold off
    text(4*size(q,2),0.47*size(q,1),'=','color','yellow','fontname','symbol','fontsize',24,'horizontalalignment','center')
    line(size(q,2)*[1 1]*0.2, [0.9 0.8]*size(q,1),'color','yellow','linewidth',1)
    line([0.2 0.4]*size(q,2),0.9*size(q,1)*[1 1],'color','yellow','linewidth',1)
    text(0.46*size(q,2),0.9*size(q,1),'k_x','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic')
    text(0.2*size(q,2),0.73*size(q,1),'k_z','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic','HorizontalAlignment','center')
    caxis([-4 0])
    h = colorbar;
    set(h,'LineWidth',0.5)
    set(h,'Fontname','calibri','FontSize',7)
    set(h,'PlotBoxAspectRatio',[1 10 1])

    print -r600 -dpng OTF4Pi
    
    % 4Pi A
    ff = zeros(size(q));
    ff(q.^2+w.^2>0.98 & q.^2+w.^2<1.02 & abs(atan(q./w))<asin(1.2/1.33) & w<0)=1;
    ff = ff*(btrans.*repmat(qv,[1 nmax]));
    ff = mConv(ff,flipud(ff));
    amp0 = ff*(btrans.*repmat(qq,[numel(qv) 1]))';
    amp0 = amp0./max(amp0(:));
    amp0(amp0<1e-6) = 0;
    ff = amp0*(btrans.*repmat(qv,[1 nmax]));
    gg = amp1*(btrans.*repmat(qv,[1 nmax]));
    ff = mConv(gg,flipud(ff));
    amp2 = ff*(btrans.*repmat(qq,[numel(qv) 1]))';
    amp2 = amp2./max(amp2(:));
    amp2(amp2<1e-6) = 0;
    
    CombineImages(log10(cat(3,[fliplr(amp1) amp1],[fliplr(amp0) amp0],[fliplr(amp2) amp2])),1,3,[],{{'(A)',''} ,{'(B)',''},{'(C)',''}})
    h = get(gca,'children');
    for j=1:3
        set(h(j),'fontname','calibri','fontsize',8)
    end
    line(2*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    line(4*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    hold on
    plot(2*size(q,2),size(q,1)*0.5,'xy','markersize',15,'linewidth',1); 
    plot(2*size(q,2),size(q,1)*0.5,'oy','markersize',15,'linewidth',1); 
    hold off
    text(4*size(q,2),0.47*size(q,1),'=','color','yellow','fontname','symbol','fontsize',24,'horizontalalignment','center')
    line(size(q,2)*[1 1]*0.2, [0.9 0.8]*size(q,1),'color','yellow','linewidth',1)
    line([0.2 0.4]*size(q,2),0.9*size(q,1)*[1 1],'color','yellow','linewidth',1)
    text(0.46*size(q,2),0.9*size(q,1),'k_x','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic')
    text(0.2*size(q,2),0.73*size(q,1),'k_z','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic','HorizontalAlignment','center')
    caxis([-4 0])
    h = colorbar;
    set(h,'LineWidth',0.5)
    set(h,'Fontname','calibri','FontSize',7)
    set(h,'PlotBoxAspectRatio',[1 10 1])
    
    print -r600 -dpng OTF4PiA
    

    
    % 4Pi C
    ff = amp1*(btrans.*repmat(qv,[1 nmax]));
    ff = mConv(ff,flipud(ff));
    amp2 = ff*(btrans.*repmat(qq,[numel(qv) 1]))';
    amp2 = amp2./max(amp2(:));
    amp2(amp2<1e-6) = 0;
    
    CombineImages(log10(cat(3,[fliplr(amp1) amp1],[fliplr(amp1) amp1],[fliplr(amp2) amp2])),1,3,[],{{'(A)',''} ,{'(B)',''},{'(C)',''}})
    h = get(gca,'children');
    for j=1:3
        set(h(j),'fontname','calibri','fontsize',8)
    end
    line(2*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    line(4*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    hold on
    plot(2*size(q,2),size(q,1)*0.5,'xy','markersize',15,'linewidth',1); 
    plot(2*size(q,2),size(q,1)*0.5,'oy','markersize',15,'linewidth',1); 
    hold off
    text(4*size(q,2),0.47*size(q,1),'=','color','yellow','fontname','symbol','fontsize',24,'horizontalalignment','center')
    line(size(q,2)*[1 1]*0.2, [0.9 0.8]*size(q,1),'color','yellow','linewidth',1)
    line([0.2 0.4]*size(q,2),0.9*size(q,1)*[1 1],'color','yellow','linewidth',1)
    text(0.46*size(q,2),0.9*size(q,1),'k_x','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic')
    text(0.2*size(q,2),0.73*size(q,1),'k_z','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic','HorizontalAlignment','center')
    caxis([-4 0])
    h = colorbar;
    set(h,'LineWidth',0.5)
    set(h,'Fontname','calibri','FontSize',7)
    set(h,'PlotBoxAspectRatio',[1 10 1])
    
    print -r600 -dpng OTF4PiC
    
end

if 1 % SIM
    
    [q,w] = meshgrid(0.0025:0.005:4,-4:0.005:4);
    amp0 = zeros(size(q));
    amp0(q.^2+w.^2>0.98 & q.^2+w.^2<1.02 & w<0 & abs(atan(q./w))<asin(1.2/1.33))=1;
    
    qq = 0:0.5:2e3; nmax = numel(qq);
    qv = q(1,:)';
    btrans = zeros(numel(qv),nmax);
    for j=1:nmax
        btrans(:,j) = besselj(0,qv/max(qv)*qq(j));
    end
    
    ff = amp0*(btrans.*repmat(qv,[1 nmax]));
    ff = mConv(ff,flipud(ff));
    amp1 = ff*(btrans.*repmat(qq,[numel(qv) 1]))';
    
    amp1 = amp1./max(amp1(:));
    amp1(amp1<1e-6) = 0;

    amp0 = zeros(size(q));
    amp0(q.^2+w.^2>0.98 & q.^2+w.^2<1.02 & w<0 & abs(atan(q./w))<asin(1.2/1.33) & abs(atan(q./w))>asin(1.1/1.33))=1;
    amp0 = [fliplr(amp0) amp0];
    amp2 = mConv2(amp0,amp0);

    amp2 = amp2./max(amp2(:));
    amp2(amp2<1e-6) = 0;
       
    CombineImages(log10(cat(3,amp0,amp0,amp2)),1,3,[],{{'(A)',''} ,{'(B)',''},{'(C)',''}})
    h = get(gca,'children');
    for j=1:3
        set(h(j),'fontname','calibri','fontsize',8)
    end
    hold
    plot(size(q,2), 0.5*size(q,1),'x','Markersize',7,'color','yellow','linewidth',1)
    hold
    line(2*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    line(4*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    hold on
    plot(2*size(q,2),size(q,1)*0.5,'xy','markersize',15,'linewidth',1); 
    plot(2*size(q,2),size(q,1)*0.5,'oy','markersize',15,'linewidth',1); 
    hold off
    text(4*size(q,2),0.47*size(q,1),'=','color','yellow','fontname','symbol','fontsize',24,'horizontalalignment','center')
    line(size(q,2)*[1 1]*0.2, [0.9 0.8]*size(q,1),'color','yellow','linewidth',1)
    line([0.2 0.4]*size(q,2),0.9*size(q,1)*[1 1],'color','yellow','linewidth',1)
    text(0.46*size(q,2),0.9*size(q,1),'k_x','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic')
    text(0.2*size(q,2),0.73*size(q,1),'k_z','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic','HorizontalAlignment','center')
    caxis([-4 0])
    h = colorbar;
    set(h,'LineWidth',0.5)
    set(h,'Fontname','calibri','FontSize',7)
    set(h,'PlotBoxAspectRatio',[1 10 1])
    
    print -r600 -dpng OTFSIMExc

    amp3 = mConv2(amp2,[fliplr(amp1) amp1]);
    amp3 = amp3./max(amp3(:));
    amp3(amp3<1e-6) = 0;

    CombineImages(log10(cat(3,amp2,[fliplr(amp1) amp1],amp3)),1,3,[],{{'(A)',''} ,{'(B)',''},{'(C)',''}})
    h = get(gca,'children');
    for j=1:3
        set(h(j),'fontname','calibri','fontsize',8)
    end
    line(2*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    line(4*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    hold on
    plot(2*size(q,2),size(q,1)*0.5,'xy','markersize',15,'linewidth',1); 
    plot(2*size(q,2),size(q,1)*0.5,'oy','markersize',15,'linewidth',1); 
    hold off
    text(4*size(q,2),0.47*size(q,1),'=','color','yellow','fontname','symbol','fontsize',24,'horizontalalignment','center')
    line(size(q,2)*[1 1]*0.2, [0.9 0.8]*size(q,1),'color','yellow','linewidth',1)
    line([0.2 0.4]*size(q,2),0.9*size(q,1)*[1 1],'color','yellow','linewidth',1)
    text(0.46*size(q,2),0.9*size(q,1),'k_x','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic')
    text(0.2*size(q,2),0.73*size(q,1),'k_z','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic','HorizontalAlignment','center')
    caxis([-4 0])
    h = colorbar;
    set(h,'LineWidth',0.5)
    set(h,'Fontname','calibri','FontSize',7)
    set(h,'PlotBoxAspectRatio',[1 10 1])

    print -r600 -dpng OTFSIM
    
end

if 1 % SPIM
    
    [q,w] = meshgrid(0.0025:0.005:4,-4:0.005:4);
    amp0 = zeros(size(q));
    amp0(q.^2+w.^2>0.98 & q.^2+w.^2<1.02 & w<0 & abs(atan(q./w))<asin(1.2/1.33))=1;
    
    qq = 0:0.5:2e3; nmax = numel(qq);
    qv = q(1,:)';
    btrans = zeros(numel(qv),nmax);
    for j=1:nmax
        btrans(:,j) = besselj(0,qv/max(qv)*qq(j));
    end
    
    ff = amp0*(btrans.*repmat(qv,[1 nmax]));
    ff = mConv(ff,flipud(ff));
    amp1 = ff*(btrans.*repmat(qq,[numel(qv) 1]))';
    
    amp1 = amp1./max(amp1(:));
    amp1(amp1<1e-6) = 0;

    amp0 = zeros(size(q));
    amp0(q.^2+w.^2>0.98 & q.^2+w.^2<1.02 & abs(atan(w./q))<asin(0.3) & abs(atan(q./w))>asin(1.1/1.33))=1;
    amp0 = [zeros(size(q)) amp0];
    amp2 = mConv2(amp0,amp0);

    amp2 = amp2./max(amp2(:));
    amp2(amp2<1e-6) = 0;
       
    CombineImages(log10(cat(3,amp0,amp0,amp2)),1,3,[],{{'(A)',''} ,{'(B)',''},{'(C)',''}})
    h = get(gca,'children');
    for j=1:3
        set(h(j),'fontname','calibri','fontsize',8)
    end
    hold
    plot(size(q,2), 0.5*size(q,1),'x','Markersize',7,'color','yellow','linewidth',1)
    hold
    line(2*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    line(4*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    hold on
    plot(2*size(q,2),size(q,1)*0.5,'xy','markersize',15,'linewidth',1); 
    plot(2*size(q,2),size(q,1)*0.5,'oy','markersize',15,'linewidth',1); 
    hold off
    text(4*size(q,2),0.47*size(q,1),'=','color','yellow','fontname','symbol','fontsize',24,'horizontalalignment','center')
    line(size(q,2)*[1 1]*0.2, [0.9 0.8]*size(q,1),'color','yellow','linewidth',1)
    line([0.2 0.4]*size(q,2),0.9*size(q,1)*[1 1],'color','yellow','linewidth',1)
    text(0.46*size(q,2),0.9*size(q,1),'k_x','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic')
    text(0.2*size(q,2),0.73*size(q,1),'k_z','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic','HorizontalAlignment','center')
    caxis([-4 0])
    h = colorbar;
    set(h,'LineWidth',0.5)
    set(h,'Fontname','calibri','FontSize',7)
    set(h,'PlotBoxAspectRatio',[1 10 1])
    
    print -r600 -dpng OTFSPIMExc

    amp3 = mConv2(amp2,[fliplr(amp1) amp1]);
    amp3 = amp3./max(amp3(:));
    amp3(amp3<1e-6) = 0;

    CombineImages(log10(cat(3,amp2,[fliplr(amp1) amp1],amp3)),1,3,[],{{'(A)',''} ,{'(B)',''},{'(C)',''}})
    h = get(gca,'children');
    for j=1:3
        set(h(j),'fontname','calibri','fontsize',8)
    end
    line(2*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    line(4*size(q,2)*[1 1],size(q,1)*[0.4 0.6],'color','k','linewidth',5)
    hold on
    plot(2*size(q,2),size(q,1)*0.5,'xy','markersize',15,'linewidth',1); 
    plot(2*size(q,2),size(q,1)*0.5,'oy','markersize',15,'linewidth',1); 
    hold off
    text(4*size(q,2),0.47*size(q,1),'=','color','yellow','fontname','symbol','fontsize',24,'horizontalalignment','center')
    line(size(q,2)*[1 1]*0.2, [0.9 0.8]*size(q,1),'color','yellow','linewidth',1)
    line([0.2 0.4]*size(q,2),0.9*size(q,1)*[1 1],'color','yellow','linewidth',1)
    text(0.46*size(q,2),0.9*size(q,1),'k_x','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic')
    text(0.2*size(q,2),0.73*size(q,1),'k_z','color','yellow','fontname','calibri','fontsize',7,'fontangle','italic','HorizontalAlignment','center')
    caxis([-4 0])
    h = colorbar;
    set(h,'LineWidth',0.5)
    set(h,'Fontname','calibri','FontSize',7)
    set(h,'PlotBoxAspectRatio',[1 10 1])

    print -r600 -dpng OTFSPIM
    
end

return % test of Bessel transform

q = 0.0005:0.001:2.5;
amp0 = zeros(size(q)); 
amp0(q.^2>0.98 & q.^2<1.02)=1;

plot(log(amp0))

nmax = 1000;
qq = besselzero(0,nmax,1)';
btrans = zeros(numel(q),nmax);
for j=1:nmax
    btrans(:,j) = besselj(0,q'/max(q)*qq(j));
end

ff = amp0*(btrans.*repmat(q',[1 nmax]));
ff = ff.^2;
amp1 = ff*(btrans.*repmat(qq,[numel(q) 1]))';

amp1 = amp1./max(amp1(:));
amp1(amp1<1e-10) = 0;
plot(log10(amp1))

[q2,w2] = meshgrid(-2.5:0.001:2.5,-2.5:0.001:2.5);
amp0 = zeros(size(q2)); 
amp0(q2.^2+w2.^2>0.98 & q2.^2+w2.^2<1.02)=1;
tst = mConv2(amp0,amp0);
tst = (tst-min(tst(:)))/(max(tst(:))-min(tst(:)));

plot(q,log10(amp1),q2((end+1)/2,:),log10(tst((end+1)/2,:)))
