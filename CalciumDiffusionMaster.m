% CalciumDiffusionMaster

NatConst

L = 10;
R = 50;
DiffConstant = 1e-5;
InfluxPoint = 0.5;
ChargeV = 0.5:0.5:4;

rho0v = [0.15 0.25 0.3];
IonStrengthv = [0.1 0.15 0.2]; 

fsize = 10;

for jrho0 = 1:length(rho0v)
    rho0 = rho0v(jrho0);
    for jIonStrength = 1:length(IonStrengthv)
        IonStrength = IonStrengthv(jIonStrength);
        Charge = 0;
        [dist0, flux0, rhov, zv, delta] = CalciumDiffusion(L, R, rho0, DiffConstant, Charge, IonStrength, InfluxPoint);
        for j=1:length(ChargeV)
            Charge = ChargeV(j);
            [dist, flux, rhov, zv, delta] = CalciumDiffusion(L, R, rho0, DiffConstant, Charge, IonStrength, InfluxPoint);

            ind = rhov<=5.5;
            tmp = log10([flipud(dist0(ind,1:end-1)); ones(1,length(zv)-1); dist(ind,1:end-1)]/AvogadroConstant*1e3);
            tmp((end+1)/2,:) = max(tmp(:));
            pcolor(zv(1:end-1),[-fliplr(rhov(ind)) 0 rhov(ind)],tmp)
            colormap jet
            shading interp;
            axis image
            axis([0 9.9 -5 5])
            patch([0 10 10 0],1.5*[-rho0 -rho0 rho0 rho0],'w')
            set(gca,'yticklabel',num2str(abs(str2num(get(gca,'yticklabel')))))
            xlabel('\itz\rm [nm]');
            ylabel('\rho [nm]');
            h = colorbar;
            set(get(h,'ylabel'),'String','lg \itc\rm [M]')
            eval(['print -dpng -r300 garp_conc_C' mint2str(j,3) '_r' mint2str(100*rho0,3) '_I' mint2str(100*IonStrength,3)])
            set(h,'yscale','log')
            eval(['print -dpng -r300 log_conc_C' mint2str(j,3) '_r' mint2str(100*rho0,3) '_I' mint2str(100*IonStrength,3)])
            
            tmp = log10(-[flipud(diff(dist0(ind,:),[],2)); -ones(1,length(zv)-1); diff(dist(ind,:),[],2)]/AvogadroConstant/diff(zv(1:2))*1e7*DiffConstant);
            tmp((end+1)/2,:) = max(tmp(:));
            pcolor(zv(1:end-1),[-fliplr(rhov(ind)) 0 rhov(ind)],real(tmp))
            shading interp
            axis image
            axis([0 9.9 -5 5])
            patch([0 10 10 0],1.5*[-rho0 -rho0 rho0 rho0],'w')
            set(gca,'yticklabel',num2str(abs(str2num(get(gca,'yticklabel')))))
            xlabel('\itz\rm [nm]');
            ylabel('\rho [nm]');
            h = colorbar;
            set(get(h,'ylabel'),'String','lg -\itD\rm\cdot\partial\itc\rm/\partial\itz\rm [mol/cm^2/s]')
            eval(['print -dpng -r300 garp_flux_C' mint2str(j,3) '_r' mint2str(100*rho0,3) '_I' mint2str(100*IonStrength,3)])
            set(h,'yscale','log')
            eval(['print -dpng -r300 log_flux_C' mint2str(j,3) '_r' mint2str(100*rho0,3) '_I' mint2str(100*IonStrength,3)])
            
            semilogy(rhov,flux/AvogadroConstant,rhov,flux0/AvogadroConstant);
            xlabel('\rho [nm]')
            ylabel('-\itD\rm\cdot\partial\itc\rm/\partial\itz\rm [mol/cm^2/s]')
            eval(['print -dpng -r300 garp_endflux_C' mint2str(j,3) '_r' mint2str(100*rho0,3) '_I' mint2str(100*IonStrength,3)])
            
        end
        %eval(['print -dpng -r300 garp_r' mint2str(100*rho0,3) '_I' mint2str(100*IonStrength,3)])
    end
end

return

surf(zv(1:end-1),[-fliplr(rhov) 0 rhov],tmp);
shading interp
set(gca,'plotboxaspectratio',[1 2*R/L 3])
ax = axis;
axis([0 L -R R -30 -5])
caxis([-30 -5])
camlight
set(gca,'gridlinestyle','--','linewidth',0.5)
view([60 20])
alpha(0.8)
text(30,-3,-30,'\rho [nm]','horizontalalignment','center','fontsize',fsize); % xlabel('\itz\rm [nm]');
text(5,-65,-30,'\itz\rm [nm]','horizontalalignment','center','fontsize',fsize); % ylabel('\rho [nm]');
zlabel('lg \itc\rm [M]','fontsize',fsize)
title([mnum2str(Charge,1,1) ' \ite\rm^-/nm'],'fontsize',1.2*fsize)
