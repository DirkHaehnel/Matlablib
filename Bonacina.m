% NAv= 0.3:0.02:0.5; for joerg=1:length(NAv) NA=NAv(joerg); Bonacina;
% title(['NA = ' mnum2str(NA,1,2)]); eval(['print -dpng -r300 BonacinaNA'
% mint2str(1e2*NA,3)]); end

close all

if 1
    %close all
    clear imtotal final
    alv = -35/180*pi;
    focposv = linspace(684,754,8);
    linebar = 10;
    vv = [1 1 1 1 1 1 1 1]; % brightness scaling
    
    for jf = 1:length(focposv)
        focpos = focposv(jf);
        maxfield = 20;
        shift = [-1.29 -6.84];
        al = pi/2;
        be = [pi/2 pi/2];
        brightness = [1 1.2];
        
        rho = [0 sqrt(maxfield^2+(maxfield+sqrt(sum(shift.^2))/2)^2)];
        z = 0;
        NA = 0.6;
        n1 = 1;
        n = 1;
        n2 = [2.092 1.52];
        %n2 = [2.092 1.52 1];
        d1 = [];
        d = 0;
        d2 = 0.159;
        %d2 = [0.159 700];
        lambda = 0.39;
        mag = 35;
        
        atf = [1.51 -2e3];
        
        [intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole(rho, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf);
        
        pixel = 7.4;
        nn = ceil(mag*[maxfield (maxfield+sqrt(sum(shift.^2))/2)]/pixel);

        [int1, x, y, ex1, ey1, bx1, by1] = SEPImage(al(1),be(1),nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
        ex1 = brightness(1)*interp2(x,y,ex1,x+mag/pixel*shift(1)/2,y+mag/pixel*shift(2)/2,'cubic',0);
        ey1 = brightness(1)*interp2(x,y,ey1,x+mag/pixel*shift(1)/2,y+mag/pixel*shift(2)/2,'cubic',0);
        bx1 = brightness(1)*interp2(x,y,bx1,x+mag/pixel*shift(1)/2,y+mag/pixel*shift(2)/2,'cubic',0);
        by1 = brightness(1)*interp2(x,y,by1,x+mag/pixel*shift(1)/2,y+mag/pixel*shift(2)/2,'cubic',0);
        for ja=1:length(alv)
            [int2, x, y, ex2, ey2, bx2, by2] = SEPImage(alv(ja),be(2),nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
            ex2 = brightness(2)*interp2(x,y,ex2,x-mag/pixel*shift(1)/2,y-mag/pixel*shift(2)/2,'cubic',0);
            ey2 = brightness(2)*interp2(x,y,ey2,x-mag/pixel*shift(1)/2,y-mag/pixel*shift(2)/2,'cubic',0);
            bx2 = brightness(2)*interp2(x,y,bx2,x-mag/pixel*shift(1)/2,y-mag/pixel*shift(2)/2,'cubic',0);
            by2 = brightness(2)*interp2(x,y,by2,x-mag/pixel*shift(1)/2,y-mag/pixel*shift(2)/2,'cubic',0);
            
            int = real((ex1*cos(be(1))-ey1*sin(be(1))+ex2*cos(be(2))-ey2*sin(be(2))).*conj(by1*cos(be(1))+bx1*sin(be(1))+by2*cos(be(2))+bx2*sin(be(2))) - ...
                (ey1*cos(be(1))+ex1*sin(be(1))+ey2*cos(be(2))+ex2*sin(be(2))).*conj(bx1*cos(be(1))-by1*sin(be(1))+bx2*cos(be(2))-by2*sin(be(2))));
            
            int1 = real((ex1*cos(be(1))-ey1*sin(be(1))).*conj(by1*cos(be(1))+bx1*sin(be(1))) - ...
                (ey1*cos(be(1))+ex1*sin(be(1))).*conj(bx1*cos(be(1))-by1*sin(be(1)))) + ...
                real((ex2*cos(be(2))-ey2*sin(be(2))).*conj(by2*cos(be(2))+bx2*sin(be(2))) - ...
                (ey2*cos(be(2))+ex2*sin(be(2))).*conj(bx2*cos(be(2))-by2*sin(be(2))));
            
            
            rr = nn(2)-nn(1);
            final(:,:,jf) = int;
            final1(:,:,jf) = int1;
            %imtotal((jf-1)*size(int,1)+1:jf*size(int,1),(ja-1)*size(int,2)+1:ja*size(int,2)) = int/max(int(:));
            imtotal(floor((jf-1)/4)*size(int,1)+1:(floor((jf-1)/4)+1)*size(int,1),rem(jf-1,4)*size(int,2)+1:(rem(jf-1,4)+1)*size(int,2)) = int/max(int(:))/vv(jf);
            imtotal1(floor((jf-1)/4)*size(int,1)+1:(floor((jf-1)/4)+1)*size(int,1),rem(jf-1,4)*size(int,2)+1:(rem(jf-1,4)+1)*size(int,2)) = int1/max(int1(:))/vv(jf);
        end
        mim(imtotal)
    end
    if ~isempty(linebar)
        mim(mConv2(imtotal,Disk(2)))
        caxis([0 0.4]);
        hold on
        ax = axis;
        line([ax(2)-15-mag*linebar/pixel ax(2)-15],[15 15],'Color','y','linewidth',3)
        hold off

        figure
        mim(mConv2(imtotal1,Disk(2)))
        caxis([0 0.4]);
        hold on
        ax = axis;
        line([ax(2)-15-mag*linebar/pixel ax(2)-15],[15 15],'Color','y','linewidth',3)
        hold off
    end
end

% return

% make experimental figure
% load C:\Joerg\MATLAB\BonacinaData.mat
figure
mim(pict(:,:,1))
[a,b] = PickPixel;
tmp = []; 
sx = (size(int,1)-1)/2; sx = -sx:sx;
sy = (size(int,2)-1)/2; sy = -sy:sy;
cnt = 1; for k=14:-1:11 tmp(1:length(sx),(cnt-1)*length(sy)+1:cnt*length(sy)) = pict(mean(a)+sx,mean(b)+sy,k)/max(max(max(pict(mean(a)+sx,mean(b)+sy,k)))); cnt=cnt+1; end;
cnt = 1; for k=10:-1:7 tmp(length(sx)+1:2*length(sx),(cnt-1)*length(sy)+1:cnt*length(sy)) = pict(mean(a)+sx,mean(b)+sy,k)/max(max(max(pict(mean(a)+sx,mean(b)+sy,k)))); cnt=cnt+1; end;
mim(tmp)
linebar = 10;
hold on
ax = axis;
line([ax(2)-15-linebar*mag/pixel ax(2)-15],[15 15],'Color','y','linewidth',3)
hold off

