bckgrnd = 0;
signalv = [5e2 1e3 2e3 4e3 8e3 16e3];
max_cnt = 1e4;

n = 1.333;
n1 = 1.333;
d0 = [];
d = 0; d1 = [];
lamem = 0.575;
mag = 250;
focus = 0;

nn = 15;
n2 = 2*nn+1;
n3 = (2*nn+1)^2;
[x,y] = meshgrid(-nn:nn,-nn:nn);
al_res = 10;
alv = (0:al_res:90)/180*pi;
coord = (-nn:nn)';

pixelv = [12.5 25 50 100];
NAv = [1.2 1.4 1.45 1.65];

for jNA=1:length(NAv)
    NA = NAv(jNA);
    n0 = 1.52;
    if NA==1.65
        n0 = 1.78;
    end
    if NA==1.2
        n0 = n;
    end
    [intx inty intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 8.75], 0, NA, n0, n, n1, d0, d, d1, lamem, mag, focus);
    rho0 = rho;
    for jpix = 1:length(pixelv)
        rho = rho0/pixelv(jpix);
        for k=1:length(alv)
            al = alv(k);
            for s=1:length(signalv)
                signal = signalv(s);
                for j=1:max_cnt
                    x0 = rand-0.5; y0 = rand-0.5; be = rand*pi*2;
                    r = sqrt((x-x0).^2+(y-y0).^2);
                    p = angle((x-x0)+1i*(y-y0));
                    int = real((cos(al)*(interp1(rho,fxx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,fxx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,fxz,r,'cubic')).*...
                        conj(cos(al)*(interp1(rho,byx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,byx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,byz,r,'cubic')) + ...
                        (cos(al)*sin(2*(p-be)).*interp1(rho,fxx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,fxz,r,'cubic')).*...
                        conj(cos(al)*sin(2*(p-be)).*interp1(rho,byx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,byz,r,'cubic')));
                    int = int/sum(int(:));
                    z = poissrnd(bckgrnd+signal*int);
                    z(isnan(z)) = 0;
                    tmp = fit(coord,sum(z)','gauss1');
                    px(j,k,s) = tmp.b1-x0;
                    tmp = fit(coord,sum(z')','gauss1');
                    py(j,k,s) = tmp.b1-y0;
                    if mod(j,100)==0 j, end
                end
            end
            eval(['save PositionAccuracy' mint2str(10*pixelv(jpix),4) 'pix' mint2str(100*NA,3) 'NA']);
        end
    end
end

return

pixelv = [12.5 25 50 100];
NAv = [1.2 1.4 1.45 1.65];
al_res = 10;
alv = (0:al_res:90)/180*pi;
signalv = [5e2 1e3 2e3 4e3 8e3 16e3];
mag = 250;
bin = [-100:0.1:100]';
clear hx hy hh
for jNA=1:length(NAv)
    NA = NAv(jNA);
    for jpix = 1:length(pixelv)
        pixel = pixelv(jpix);
        eval(['load C:\Joerg\Matlab\PositionAccuracy' mint2str(10*pixelv(jpix),4) 'pix' mint2str(100*NA,3) 'NA px py']);
        for k=1:size(px,2)
            for s=1:size(px,3)
                hx(:,k,s,jpix,jNA)=mHist(pixel*px(:,k,s)/mag*1e3,bin);
                hy(:,k,s,jpix,jNA)=mHist(pixel*py(:,k,s)/mag*1e3,bin);
                hh(:,k,s,jpix,jNA) = hx(:,k,s,jpix,jNA)+hy(:,k,s,jpix,jNA);
                mm(k,s,jpix,jNA) = sum(hh(:,k,s,jpix,jNA).*bin)/sum(hh(:,k,s,jpix,jNA));
                sig(k,s,jpix,jNA) = sum(hh(:,k,s,jpix,jNA).*bin.^2)/sum(hh(:,k,s,jpix,jNA)) - mm(k,s,jpix,jNA).^2;
            end
        end
    end
end
save PositionAccuracy NAv pixelv alv signalv hx hy hh mm sig bin mag 

return

% load PositionAccuracy
% bin = [-100:0.1:100]';
% clear z
% for jNA=1:length(NAv)
%     for jpix = 1:length(pixelv)
%         for k=1:length(alv)
%             for s=1:length(signalv)
%                 ind = hh(:,k,s,jpix,jNA)>0;
%                 tst = bin(hh(:,k,s,jpix,jNA)==max(hh(:,k,s,jpix,jNA)));
%                 p(:,k,s,jpix,jNA) = simplex('GaussDouble',abs([tst(1) 0 sqrt(sig(k,s,jpix,jNA)) sqrt(sig(k,s,jpix,jNA))]),...
%                     [],[],[],[],bin(ind),hh(ind,k,s,jpix,jNA)/sum(hh(:,k,s,jpix,jNA))); 
%                 [bla, bla, z(:,k,s,jpix,jNA)] = GaussDouble(p(:,k,s,jpix,jNA),bin,hh(:,k,s,jpix,jNA)/sum(hh(:,k,s,jpix,jNA)));
%             end
%         end
%     end
% end


return

load PositionAccuracy
sig = sqrt(sig);
for jNA=1:1%length(NAv)
    NA = NAv(jNA);
    for jpix = 1:1%length(pixelv)
        pixel = pixelv(jpix)*1e3/mag;
        plot(0:90,[interp1(0:10:90,sig(:,1,jpix,jNA),0:90,'cubic')' ...
            interp1(0:10:90,sig(:,2,jpix,jNA),0:90,'cubic')' ...
            interp1(0:10:90,sig(:,3,jpix,jNA),0:90,'cubic')' ...
            interp1(0:10:90,sig(:,4,jpix,jNA),0:90,'cubic')' ...
            interp1(0:10:90,sig(:,5,jpix,jNA),0:90,'cubic')' ...
            interp1(0:10:90,sig(:,6,jpix,jNA),0:90,'cubic')']);
        xlabel('inclination angle [°]')
        ylabel('\sigma [nm]')
        title(['N.A. = ' num2str(NA) ' & pixel size = ' num2str(pixel) ' nm'])
        eval(['print -dpng PositionAccuracy' mint2str(1e3/mag*pixelv(jpix),3) 'pix' mint2str(100*NA,3) 'NA']);
    end
end


load PositionAccuracy
sig = sqrt(sig);
set(0,'defaultAxesFontName', 'Times', 'defaultAxesFontSize', 6, 'defaultLineLineWidth', 1, ...
    'defaultTextFontName', 'Times', 'defaultTextFontSize', 6);
set(0,'defaultFigurePosition',[20 34 909 913])
set(0,'defaultaxeslinewidth',1)
set(0,'defaultLineMarkerSize',1)
for jNA=1:length(NAv)
    NA = NAv(jNA);
    for jpix = 1:3%length(pixelv)
        subplot(4,4,jpix+(jNA-1)*4)
        pixel = pixelv(jpix)*1e3/mag;
        plot(90:-10:0,sig(:,:,jpix,jNA),'-o');
        colorize(get(gca,'children'));
        xlabel('inclination angle [°]')
        ylabel('\sigma [nm]')
        title(['N.A. = ' num2str(NA) ' & pixel size = ' num2str(pixel) ' nm'])
        ax = axis;
        axis([0 90 ax(3) ax(4)])
        % eval(['print -dpng -r600 PositionAccuracy' mint2str(1e3/mag*pixelv(jpix),3) 'pix' mint2str(100*NA,3) 'NA']);
    end
end


return

% Dickson's problem

NA = 1.4;
n0 = 1.52;
n = 1.333;
n1 = 1.333;
d0 = [];
d = 0; d1 = [];
lamem = 0.575;
mag = 250;
focus = 0;

nn = 15;
n2 = 2*nn+1;
n3 = (2*nn+1)^2;
[x,y] = meshgrid(-nn:nn,-nn:nn);
al_res = 10;
alv = (0:al_res:90)/180*pi;
coord = (-nn:nn)';
pixel = 2.5*4.5;
NAv = 1.4;

[intx inty intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 3], 0, NA, n0, n, n1, d0, d, d1, lamem, mag, focus);
rho = rho/pixel;

clear int
for k=1:length(alv)
    al = alv(k);
    x0 = 0.; y0 = 0.; be = 0;
    r = sqrt((x-x0).^2+(y-y0).^2);
    p = angle((x-x0)+i*(y-y0));
    int(:,:,k) = real((cos(al)*(interp1(rho,fxx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,fxx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,fxz,r,'cubic')).*...
        conj(cos(al)*(interp1(rho,byx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,byx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,byz,r,'cubic')) + ...
        (cos(al)*sin(2*(p-be)).*interp1(rho,fxx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,fxz,r,'cubic')).*...
        conj(cos(al)*sin(2*(p-be)).*interp1(rho,byx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,byz,r,'cubic')));
    subplot(2,5,k); pcolor(int(:,:,k)); shading interp; axis image; axis off; 
end
figure
plot((-nn:nn)*pixel*1e3/mag,squeeze(sum(int)))

f = fitoptions('gauss1','startpoint',[3000 0 100]);
for k=1:length(alv)
    tmp = fit((-nn:nn)'*pixel*1e3/mag,16000*squeeze(sum(int(:,:,k))/sum(sum(int(:,:,k))))','gauss1',f);
    subplot(5,2,mod(2*k-1,10)+floor((2*k-1)/10)); 
    plot((-nn:nn)'*pixel*1e3/mag,16000*squeeze(sum(int(:,:,k))/sum(sum(int(:,:,k))))',...
        (-nn:nn)'*pixel*1e3/mag,tmp.a1*exp(-(((-nn:nn)'*pixel*1e3/mag-tmp.b1)/tmp.c1).^2))
    hold on
	plot([1 1]*tmp.b1,[0 3000],':','linewidth',1)
    hold off
    set(gca,'fontsize',10) 
    axis([-600 600 0 3000])
    text(300,2500,[mnum2str(tmp.b1,2,1) ' nm'],'fontsize',12)
    text(-550,2500,[int2str(10*(k-1)) '°'],'fontsize',12)
end

return

for k=1:length(alv)
    [mx(k), my(k)] = GaussEllipse((-nn:nn)*pixel*1e3/mag, (-nn:nn)*pixel*1e3/mag, 16000*squeeze(int(:,:,k))/sum(sum(int(:,:,k))));
end

