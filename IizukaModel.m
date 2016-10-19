function mask = IizukaModel()

close all

NA = 1.45; % numerical aperture
n0 = 1.52; % immersion medium
n = 1.38; % sample medium
n1 = n;
d0 = [];
d = 0; d1 = [];
lamex = 0.570; % emission wavelength
mag = 300; % magnification
focusv = (0.95:0.02:1.01); % defcousing in mum
%focusv = -(0.9:0.025:1); % defcousing in mum
pixel = 16; % CCD pixel size in mum

nn = 15; % half size of mask size for pattern matching

al_res = 10; % minimum resolution of out-of-plane angle
theta = (90:-al_res:0)/180*pi;

for jv=1:length(focusv)
    focus = focusv(jv);
    [intx inty intz, rho, tmp, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 ceil(nn*sqrt(2)*pixel/mag)], 0, NA, n0, n, n1, d0, d, d1, lamex, mag, focus);
    for cnt=1:length(theta)
        mask = SEPImage(theta(cnt),0,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
        subplot(length(focusv),length(theta),(jv-1)*length(theta) + cnt);
        mim(mask)
        line([0.5 size(mask,2)+0.5 size(mask,2)+0.5 0.5 0.5],[0.5 0.5 size(mask,1)+0.5 size(mask,1)+0.5 0.5],'color','k','linewidth',1)
        text(nn-7,-nn/7,['\theta = ' mint2str(theta(cnt)/pi*180) '°'],'fontsize',3,'fontname','times');
    end
    %colormap(flipud(gray))
end
%print -dpng -r300 ModelPatterns



