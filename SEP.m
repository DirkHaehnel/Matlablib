if 1
    NA = 1.49;
    focposv = 0:0.2:1;
    cnt = 1;
    for jf = 1:length(focposv)
        focpos = focposv(jf);   
        [intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 5], 0, NA, 1.51, 1.33, 1.33, [], 0, [], 0.57, 240, focpos); 
        al = pi/2; be = 0; 
        int = SEPImage(al,be,20,17,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
        subplot(2,6,cnt); 
        mim(int)        
        title([mnum2str(focpos,1,1) ' \mum'],'fontname','times','fontsize',10);
        drawnow
        be = pi/2; al = 0; 
        int = SEPImage(al,be,20,17,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
        subplot(2,6,cnt+6); 
        mim(int)        
        title([mnum2str(focpos,1,1) ' \mum'],'fontname','times','fontsize',10);
        drawnow
        cnt = cnt+1;
    end
end

if 0
    NA = 1.24;
    al = pi/2; be = 0;
    thickness = 100;
    focposv = thickness:0.5:thickness+50;
    for j=1:length(focposv)
        focpos = focposv(j);
        [intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 0], 0, NA, [1.52 1.46], 1.0, 1.0, thickness, 0, [], 0.58, 100, focpos); 
        res(j) = intx(1);
    end
    f0 = focposv(res==max(res));
    al = pi/2; be = 0; cnt = 1;
    for focpos = 0.1:0.1:2
        [int, rho, phi] = SEPCircularScatter([0 2], 0, NA, [1.52 1.46], 1.0, 1.0, thickness, 0, [], 0.58, 100, f0 + focpos); 
        subplot(4,5,cnt); 
        pcolor(cos(phi)*rho, sin(phi)*rho, int); axis image; shading interp; axis off; colormap hot        
        title([mnum2str(focpos,1,1) ' \mum'],'fontname','times','fontsize',10);
        drawnow
        cnt = cnt+1;
    end
end


if 0
    al = pi/2; be = 0; cnt = 1;
    for focpos = 0.1:0.1:1.5
        NA = 1.24;
        [intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 2], 0, NA, 1.52, 1.0, 1.0, [], 0, [], 0.58, 100, 5 + focpos, [1.46 100]); 
        inta = SEPImage(al,be,10,6.5,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
        intb = SEPImage(al,be+pi/2,10,6.5,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
        subplot(3,5,cnt); 
        mim(inta+intb)        
        title([mnum2str(focpos,1,1) ' \mum'],'fontname','times','fontsize',10);
        drawnow
        cnt = cnt+1;
    end
end

if 0
    al = 0; mm = 5;
    
    for focpos = 0.5:0.1:4.0
        for kk=1:3
            for j = 1:11
                alpha = [zeros(1,mm) 2*pi/20*(j-6)/5];
                for k = 1:11
                    NA = 1.2 + 0.1*(kk-1) + (k-1)*0.01
                    [intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 2], 0, NA, 1.52, 1.52, 1, [], 0, [], 0.58, 100, focpos, [], alpha); 
                    col = cos(0*phi);
                    int = real((cos(al)*(col*fxx0+cos(2*phi)*fxx2)+sin(al)*cos(phi)*fxz).*conj(cos(al)*(col*byx0+cos(2*phi)*byx2)+sin(al)*cos(phi)*byz) + (cos(al)*sin(2*phi)*fxx2+sin(al)*sin(phi)*fxz).*conj(cos(al)*sin(2*phi)*byx2+sin(al)*sin(phi)*byz));
                    subplot(11,11,(j-1)*11+k); 
                    pcolor(cos(phi)*rho,sin(phi)*rho,int); 
                    axis off; 
                    axis image; 
                    if j==1
                        title(mnum2str(NA,1,2),'fontname','times','fontsize',10);
                    end
                    if k==1
                        h = text(-250,0,mnum2str(alpha(end)/pi/2,1,3),'fontname','times','fontsize',10);
                        set(h,'horizontalalignment','right');
                    end
                    shading interp; 
                    colormap gray; 
                    drawnow
                end; 
            end
            eval(['print -dpng aberration' mint2str(length(alpha),2) 'p' mint2str(10*focpos,2) 'z' mint2str(10*(NA-1)-1) mint2str(10*(NA-1))]);
        end
    end
    
end