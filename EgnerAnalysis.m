if ~exist('ch1')
    % Die Reihenfolge in z ist: ch1->ch3->ch4->ch2
    load c:\Joerg\Doc\Microscopy\SOFI\Egner\SOFI_02\mitos_trafo.mat 
    ch1 = double(ch1);
    tmp = double(ch2);
    ch2 = double(ch3);
    ch3 = double(ch4);    
    ch4 = tmp;
end

[sof1, im1] = SOFIXAnalysis(ch1);
[sof2, im2] = SOFIXAnalysis(ch2);
[sof3, im3] = SOFIXAnalysis(ch3);
[sof4, im4] = SOFIXAnalysis(ch4);
[sof12, im12] = SOFIXAnalysis(squeeze(interp1(1:2,permute(cat(4,ch1,ch2),[4 1 2 3]),1.5,'spline')));
[sof23, im23] = SOFIXAnalysis(squeeze(interp1(1:2,permute(cat(4,ch2,ch3),[4 1 2 3]),1.5,'spline')));
[sof34, im34] = SOFIXAnalysis(squeeze(interp1(1:2,permute(cat(4,ch3,ch4),[4 1 2 3]),1.5,'spline')));

save c:\Joerg\Doc\Microscopy\SOFI\Egner\SOFI_02\mitos_trafo_SOFI sof1 sof2 sof3 sof4 sof12 sof23 sof34 im1 im2 im3 im4 im12 im23 im34

j=1;
stck2 = cat(3,sof1(:,:,j),sof12(:,:,j),sof2(:,:,j),sof23(:,:,j),sof3(:,:,j),sof34(:,:,j),sof4(:,:,j));
stck2 = stck2/max(stck2(:));

% Für das Fourier-apodizieren lauten die Paramter: Pixelgröße 133nm, NA=1.2
% H2O-immersion. -> Brechungsindex in der Probe ~1.334

NA = 1.2; % numerical aperture
n1 = 1.334; % ref. index of sample
n = n1; 
n2 = n1;
d1 = [];
d = 0;
d2 = [];
lambda = 0.67; % emission wavelength in micron
mag = 60; % magnification
pix = mag*0.133/4; % virtual pixel size in micron
stck2f = SOFIXFourier(stck2,NA,n1,n,n2,d1,d,d2,lambda,mag,pix);
stck2f = stck2f/max(stck2f(:));

if 0
    [x,y,z] = meshgrid(0.5*0.133*(1:size(stck2f,1)),0.5*0.133*(1:size(stck2f,2)),lambda/4*(1:size(stck2f,3)));
    lev = (1:0.1:3)*mean(stck2f(:));
    bord = 10;
    for k=1:length(lev)
        close all
        p = patch(isosurface(x(bord+1:end-bord,bord+1:end-bord,:),y(bord+1:end-bord,bord+1:end-bord,:),z(bord+1:end-bord,bord+1:end-bord,:),stck2f(bord+1:end-bord,bord+1:end-bord,:),lev(k)));
        isonormals(x(bord+1:end-bord,bord+1:end-bord,:),y(bord+1:end-bord,bord+1:end-bord,:),z(bord+1:end-bord,bord+1:end-bord,:),stck2f(bord+1:end-bord,bord+1:end-bord,:),p)
        set(p,'FaceColor','red','EdgeColor','none');
        daspect([1 1 1])
        view([0.01 0 1]);
        if k==1
            axis tight
            ax = axis;
        else
            axis(ax)
        end
        box on
        camlight
        lighting gouraud
        eval(['print -dpng -r100 test' mint2str(k,2)]);
        if k>1 && k<length(lev)
            eval(['print -dpng -r100 test' mint2str(2*length(lev)-k,2)]);
        end
    end
end

if 0
    bord = 10;
    for k=1:size(stck2f,3)
        close all
        mim(stck2f(bord+1:end-bord,bord+1:end-bord,k))
        ax = axis;
        text(0.9*ax(2),0.05*ax(4),[int2str(k*lambda/8*1e3/n1) ' nm'],'color','y','horizontalalignment','right')
        hold on
        plot(0.05*ax(2)*[1 1],[0.05*ax(4) 0.05*ax(4)+round(mag/pix)],'y','linewidth',3)
        hold off
        eval(['print -dpng -r100 test' mint2str(k,2)]);
        if k>1 && k<size(stck2f,3)
            eval(['print -dpng -r100 test' mint2str(2*size(stck2f,3)-k,2)]);
        end
    end
end
