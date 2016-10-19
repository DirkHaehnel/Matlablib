clf

alv = (0:10:90)/180*pi;
be = 0; 

xx = -10:10;
for j=1:5
    for k=1:2
        al = alv((k-1)*5+j);
        
        subplot(6,5,(k-1)*15+j);
        Pfeil([0 0 0],1.5*[cos(be)*sin(al) sin(be)*sin(al) cos(al)])
        line([0 0],[-1 1],[0 0],'color','k','linewidth',1)
        axis image
        axis([-1.4 1.4 -0.4 0.4 -0.4 1.4]);
        axis off
        colormap jet
        shading interp
        light
        view([0 90]);

        subplot(6,5,[(k-1)*15+j+5,(k-1)*15+j+10]);

        % int = DefocImage(10, 6.5, 0, 1.3, 1.52, 1., 1., [], 0, [], 0.67, 100, 0, [], [], [al be]);
        int = DefocImage(10, 6.5, 0, 1.45, 1.52, 1.33, 1.33, [], 0, [], 0.67, 100, 0, [], [], [al be]);
        mx(k,j) = Gauss2D([],[],int);
        [ind,ind] = max(sum(int)); 
        mxx(k,j) = xx(ind)*6.5;
        axis off
        hold on
        plot([11.5 11.5],[0 21],'y','linewidth',1);
        %plot([mx mx],[0 3],'w');
        hold off
    end
end
colormap hot


