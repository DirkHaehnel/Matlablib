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

for jNA=2:3;%1:length(NAv)
    NA = NAv(jNA);
    n0 = 1.52;    
    if NA==1.65
        n0 = 1.78;
    end
    if NA==1.2
        n0 = n;
    end
    n0 
    [intx inty intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 8.75], 0, NA, n0, n, n1, d0, d, d1, lamem, mag, focus);
    rho0 = rho;
    for jpix = 1:length(pixelv)
        rho = rho0/pixelv(jpix);
        for k=1:length(alv)
            al = alv(k);
            j = 1;
            for jx=0.05:0.1:1
                for jy=0.05:0.1:1
                    for jbe=0:pi/18:pi/2
                        x0 = jx-0.5; y0 = jy-0.5; be = jbe;
                        r = sqrt((x-x0).^2+(y-y0).^2);
                        p = angle((x-x0)+i*(y-y0));
                        int = real((cos(al)*(interp1(rho,fxx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,fxx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,fxz,r,'cubic')).*...
                            conj(cos(al)*(interp1(rho,byx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,byx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,byz,r,'cubic')) + ...
                            (cos(al)*sin(2*(p-be)).*interp1(rho,fxx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,fxz,r,'cubic')).*...
                            conj(cos(al)*sin(2*(p-be)).*interp1(rho,byx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,byz,r,'cubic')));
                        z = int/sum(int(:));
                        tmp = fit(coord,sum(z)','gauss1');
                        px(j,k) = tmp.b1-x0;
                        tmp = fit(coord,sum(z')','gauss1');
                        py(j,k) = tmp.b1-y0;
                        j = j+1;
                    end
                end
            end
            [jNA jpix k]
            eval(['save PositionAccuracyInf' mint2str(10*pixelv(jpix),4) 'pix' mint2str(100*NA,3) 'NA']);
        end
    end
end

return 

pixelv = [12.5 25 50 100];
NAv = [1.2 1.4 1.45 1.65];
al_res = 10;
alv = (0:al_res:90)/180*pi;
mag = 250;
bin = [-100:0.1:100]';
clear hx hy hh
for jNA=1:length(NAv)
    NA = NAv(jNA);
    for jpix = 1:length(pixelv)
        pixel = pixelv(jpix);
        eval(['load D:\Joerg\Matlab\PositionAccuracyInf' mint2str(10*pixelv(jpix),4) 'pix' mint2str(100*NA,3) 'NA px py']);
        px = [px; -px]; py = [py; -py];
        for k=1:size(px,2)
            hx(:,k,jpix,jNA) = mHist(pixel*px(:,k)/mag*1e3,bin);
            hy(:,k,jpix,jNA) = mHist(pixel*py(:,k)/mag*1e3,bin);
            hh(:,k,jpix,jNA) = hx(:,k,jpix,jNA)+hy(:,k,jpix,jNA);
            mminf(k,jpix,jNA) = sum(hh(:,k,jpix,jNA).*bin)/sum(hh(:,k,jpix,jNA));
            siginf(k,jpix,jNA) = sum(hh(:,k,jpix,jNA).*bin.^2)/sum(hh(:,k,jpix,jNA)) - mminf(k,jpix,jNA).^2;
        end
    end
end
save PositionAccuracyInf mminf siginf

return

load PositionAccuracy
sig = sqrt(sig);
load PositionAccuracyInf
siginf = sqrt(siginf);
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
        plot(90:-10:0,sig(:,:,jpix,jNA),'-o',90:-10:0,siginf(:,jpix,jNA),'-s');
        colorize(get(gca,'children'));
        xlabel('inclination angle [°]')
        ylabel('\sigma [nm]')
        title(['N.A. = ' num2str(NA) ' & pixel size = ' num2str(pixel) ' nm'])
        ax = axis;
        axis([0 90 ax(3) ax(4)])
        % eval(['print -dpng -r600 PositionAccuracy' mint2str(1e3/mag*pixelv(jpix),3) 'pix' mint2str(100*NA,3) 'NA']);
    end
end
