% EM-Wave

len = 5;

for k=1:30
    phi = 2*pi/30*(k-1);
    
    Pfeil([0 0 -0*cos(-phi)],[0 0 cos(-phi)]); 
    hold on; 
    
    Pfeil([-0*cos(-phi) 0 0],[cos(-phi) 0 0]);
    
    for j=pi/5:pi/5:2*pi
        y = cos(j-phi);
        Pfeil([0 j/pi*len -0*y],[0 j/pi*len y]); 
        Pfeil([-0*y j/pi*len 0],[y j/pi*len 0 ]); 
    end
    y = 0:pi/100:2*pi; row = zeros(size(y)); 
    plot3(row,y/pi*len,cos(y-phi),'r','linewidth',1.2)
    plot3(cos(y-phi),y/pi*len,row,'b','linewidth',1.2)
%     plot3(row,y/pi*len,cos(y-phi),'r',row,y/pi*len,-cos(y-phi),'r','linewidth',1.2)
%     plot3(cos(y-phi),y/pi*len,row,'b',-cos(y-phi),y/pi*len,row,'b','linewidth',1.2)
    view([50 20])
    caxis([-1.5 2.5]); 
    alpha(0.9)
    axis([-1 1 0 2*len -1 1])
    axis on
    set(gca,'linewidth',1,'ytick',[0 0.25 0.5 0.75 1]*2*len,'yticklabel','');
    s = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    for j=1:5
        text(1.4,(j-1)*0.5*len,-1,s{j});
    end
    
    hold off
    eval(['print -dpng EMWave' mint2str(k,2)])
end

