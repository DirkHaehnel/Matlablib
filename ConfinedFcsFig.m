if 0 % Fig3

    load D:\Doc\FCS\ConfinedFCS\ConfinedFCSInt.mat
    rintv = [0.1 0.2 0.5 1 2 3];

    load D:\Doc\FCS\ConfinedFCS\ConfinedFCSRes.mat
    jint = 3:5;
    s = {'a','b','c'};

    for j=1:3
        subplot(121)
        tmp = int(:,:,jint(j))/int0;
        ind = tmp(:,1)>0.01;
        pcolor((2+rintv(jint(j)))*rv(ind)'*sin(pv),(2+rintv(jint(j)))*rv(ind)'*cos(pv),tmp(ind,:));
        colormap hot
        hold on;
        pcolor(-(2+rintv(jint(j)))*rv(ind)'*sin(pv),(2+rintv(jint(j)))*rv(ind)'*cos(pv),tmp(ind,:));
        pcolor((2+rintv(jint(j)))*rv(ind)'*sin(pv),-(2+rintv(jint(j)))*rv(ind)'*cos(pv),tmp(ind,:));
        pcolor(-(2+rintv(jint(j)))*rv(ind)'*sin(pv),-(2+rintv(jint(j)))*rv(ind)'*cos(pv),tmp(ind,:));
        axis image;
        shading interp
        plot(rintv(jint(j))*sin(0:pi/100:2*pi),rintv(jint(j))*cos(0:pi/100:2*pi),'y'); hold off
        colorbar
        axis off
        %eval(['print -dpng -r300 ConfinedInt' mint2str(10*radiusv(j),2)]);

        subplot(122)
        tau(1,2:end,j)=tau(1,1,j);
        pcolor((2+radiusv(j))*rv(ind)'*sin(pv),(2+radiusv(j))*rv(ind)'*cos(pv),tau(ind,:,j)/tau0);
        hold on;
        pcolor(-(2+radiusv(j))*rv(ind)'*sin(pv),(2+radiusv(j))*rv(ind)'*cos(pv),tau(ind,:,j)/tau0);
        pcolor((2+radiusv(j))*rv(ind)'*sin(pv),-(2+radiusv(j))*rv(ind)'*cos(pv),tau(ind,:,j)/tau0);
        pcolor(-(2+radiusv(j))*rv(ind)'*sin(pv),-(2+radiusv(j))*rv(ind)'*cos(pv),tau(ind,:,j)/tau0);
        axis image;
        shading interp
        plot(radiusv(j)*sin(0:pi/100:2*pi),radiusv(j)*cos(0:pi/100:2*pi),'y'); hold off
        colorbar
        patch([-0.5 0.5 0.5 -0.5],[-radiusv(j)-1.55 -radiusv(j)-1.55 -radiusv(j)-1.45 -radiusv(j)-1.45],'w')
        axis off
        
        eval(['print -dtiff -r300 Fig3' s{j}]);
    end

end

if 0 % Fig3d
    load D:\Doc\FCS\ConfinedFCS\ConfinedFCSRes.mat

    j = 2;
    semilogx(autotime,res(:,1:1:31,1,j))
    colorize
    axis([1 1e6 0 1])
    xlabel('time [a.u.]');
    ylabel('autocorrelation')
    
end

if 0 % Fig4

    load D:\Doc\FCS\ConfinedFCS\ConfinedFCSInt.mat
    rintv = [0.1 0.2 0.5 1 2 3];

    load D:\Doc\FCS\ConfinedFCS\ExcludedFCSRes.mat
    jint = 1:4;
    s = {'a','b','c'};

    for j=2:4
        subplot(121)
        tmp = int(:,:,jint(j))/int0;
        ind = tmp(:,1)>0.01;
        tmp = 1-tmp;
        pcolor((2+rintv(jint(j)))*rv(ind)'*sin(pv),(2+rintv(jint(j)))*rv(ind)'*cos(pv),tmp(ind,:));
        colormap hot
        hold on;
        pcolor(-(2+rintv(jint(j)))*rv(ind)'*sin(pv),(2+rintv(jint(j)))*rv(ind)'*cos(pv),tmp(ind,:));
        pcolor((2+rintv(jint(j)))*rv(ind)'*sin(pv),-(2+rintv(jint(j)))*rv(ind)'*cos(pv),tmp(ind,:));
        pcolor(-(2+rintv(jint(j)))*rv(ind)'*sin(pv),-(2+rintv(jint(j)))*rv(ind)'*cos(pv),tmp(ind,:));
        axis image;
        shading interp
        plot(rintv(jint(j))*sin(0:pi/100:2*pi),rintv(jint(j))*cos(0:pi/100:2*pi),'y'); hold off
        colorbar
        axis off
        %eval(['print -dpng -r300 ConfinedInt' mint2str(10*radiusv(j),2)]);

        subplot(122)
        tau(1,2:end,j)=tau(1,1,j);
        pcolor((2+radiusv(j))*rv(ind)'*sin(pv),(2+radiusv(j))*rv(ind)'*cos(pv),tau(ind,:,j)/tau0);
        hold on;
        pcolor(-(2+radiusv(j))*rv(ind)'*sin(pv),(2+radiusv(j))*rv(ind)'*cos(pv),tau(ind,:,j)/tau0);
        pcolor((2+radiusv(j))*rv(ind)'*sin(pv),-(2+radiusv(j))*rv(ind)'*cos(pv),tau(ind,:,j)/tau0);
        pcolor(-(2+radiusv(j))*rv(ind)'*sin(pv),-(2+radiusv(j))*rv(ind)'*cos(pv),tau(ind,:,j)/tau0);
        axis image;
        shading interp
        plot(radiusv(j)*sin(0:pi/100:2*pi),radiusv(j)*cos(0:pi/100:2*pi),'y'); hold off
        colorbar
        patch([-0.5 0.5 0.5 -0.5],[-radiusv(j)-1.55 -radiusv(j)-1.55 -radiusv(j)-1.45 -radiusv(j)-1.45],'w')
        axis off
        
        eval(['print -dtiff -r300 Fig4' s{j-1}]);
    end

end

if 1

    clear y y0 t
    names = dir('D:\Doc\Fcs\ConfinedFCS\061027_z-scan_Atto655+GUV\*_b.mat')
    for j=1:10 [y(:,:,j) t]=FCSCrossRead(['D:\Doc\Fcs\ConfinedFCS\061027_z-scan_Atto655+GUV\' names(j).name]); end
    y0=FCSCrossRead('D:\Doc\Fcs\ConfinedFCS\061027_z-scan_Atto655+GUV\Loesung_510.mat');

    y = squeeze(sum(y(:,1:2,:),2));
    y0 = sum(y0(:,1:2),2);

    for j=1:size(y,2)
        y(:,j) = y(:,j)-mean(y(end-10:end,j));
        y(:,j) = y(:,j)./mean(y(1:10,j));
    end
    y0 = y0-mean(y0(end-10:end));
    y0 = y0./mean(y0(1:10));

    semilogx(t,y0,'oc');
    hold on
    semilogx(t,y(:,3)*1.15,t,y(:,[5 10]))
    hold off
    colorize
    ax = axis;
    axis([ax(1) max(t) 0 ax(4)])
    xlabel('time [s]')
    ylabel('autocorrelation [a.u.]')

end