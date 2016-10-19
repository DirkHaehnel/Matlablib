% Quantum Dot Images

close all
focus = 0.1:0.1:4;

mag = 110;
NA = 1.2;
lambda = 0.57;

if 0 % simple dipole
    for k=1:length(focus)
        [intx inty int, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
            SEPDipole([0 3], 0, NA, 1.51, 1., 1., [], 0, [], lambda, mag, focus(k));
        subplot(5,8,k)
        pcolor(cos(phi)*rho,sin(phi)*rho,intx+inty);axis image; shading interp; axis off
        title(num2str(focus(k)),'fontsize',12);
        colormap(gray)
        drawnow
    end
end

if 1 % circular 2D-emitter + 1D-emitter
    clear mask
    nn = 20;
    n2 = 2*nn;
    n3 = (2*nn+1)^2;
    [x,y] = meshgrid(-nn:nn,-nn:nn);
    p = angle(x+i*y);
    r = sqrt(x.^2+y.^2);
    
    for focus = 1.05:0.05:1.25;
        cnt = 1;
        for k=1:7
            om = pi/12*(k-1);
            [tmp, rho, phi, intc1, ints1] = QDSEP(0,0,om,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus);
            [tmp, rho, phi, intc2, ints2] = QDSEP(0,pi/2,om,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus,pi/2);
            rho = rho/6.45;
            int = interp1(rho,real(intc1(1,:)+intc2(1,:)),r,'cubic');
            for jj=1:4 
                int = int + cos(jj*p).*interp1(rho,real(intc1(jj+1,:)),r,'cubic') + ...
                    sin(jj*p).*interp1(rho,real(ints1(jj,:)),r,'cubic'); 
                int = int + cos(jj*(p+pi/2)).*interp1(rho,real(intc2(jj+1,:)),r,'cubic') + ...
                    sin(jj*(p+pi/2)).*interp1(rho,real(ints2(jj,:)),r,'cubic');                 
            end
            [tmp, tmp, tmp, intc1, ints1] = QDSEP(-1,0,om+pi/2,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus);
            [tmp, tmp, tmp, intc2, ints2] = QDSEP(-1,pi/2,om+pi/2,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus,pi/2);
            tmp = interp1(rho,real(intc1(1,:)+intc2(1,:)),r,'cubic');
            for jj=1:4 
                tmp = tmp + cos(jj*p).*interp1(rho,real(intc1(jj+1,:)),r,'cubic') + ...
                    sin(jj*p).*interp1(rho,real(ints1(jj,:)),r,'cubic'); 
                tmp = tmp + cos(jj*(p+pi/2)).*interp1(rho,real(intc2(jj+1,:)),r,'cubic') + ...
                    sin(jj*(p+pi/2)).*interp1(rho,real(ints2(jj,:)),r,'cubic');                 
            end
            for j=1:11
                kappa = (j-1)/10;        
                mask(:,:,cnt) = (1-kappa)*int + kappa*tmp;
                subplot(7,11,cnt)
                mim(mask(:,:,cnt));
                if k==1
                    title(mnum2str(kappa,1,2),'fontname','times','fontsize',10);
                end
                if j==1
                    h = text(-10,21,int2str(om/pi*180),'fontname','times','fontsize',10);
                    set(h,'horizontalalignment','right');
                end
                drawnow
                cnt = cnt + 1;
            end
        end
        eval(['print -dpng -r300 mask2+1d120NA' int2str(100*focus) 'foc'])
    end
end

if 0 % circular 2D-emitter + 1D-emitter
    for focus = 1:0.05:1.25
        cnt = 1;
        for k=1:7
            om = pi/12*(k-1);
            [int1, rho, phi, intc1, ints1] = QDSEP(0,0,om,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus);
            int2 = QDSEP(0,pi/2,om,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus,pi/2);
            int3 = QDSEP(-1,0,om+pi/2,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus);
            int4 = QDSEP(-1,pi/2,om+pi/2,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus,pi/2);
            for j=1:11
                kappa = (j-1)/10;        
                subplot(7,11,cnt)
                pcolor(cos(phi)*rho,sin(phi)*rho,(1-kappa)*(int1+int2) + kappa*(int3+int4));
                axis image; shading interp; axis off
                colormap hot
                if k==1
                    title(mnum2str(kappa,1,2),'fontname','times','fontsize',10);
                end
                if j==1
                    h = text(-250,11,int2str(om/pi*180),'fontname','times','fontsize',10);
                    set(h,'horizontalalignment','right');
                end
                drawnow
                cnt = cnt + 1;
            end
        end
        eval(['print -dpng -r300 q2+1d120NA' int2str(100*focus) 'foc'])
    end
end

if 0 % elliptic 2D-emitter
    for j=1:11
        kappa = -1+(j-1)/10*2;        
        for k=1:6
            om = pi/10*(k-1);
            for focus = 0.8:0.1:1.2
                [int1, rho, phi] = QDSEP(kappa,0,om,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus);
                int2 = QDSEP(kappa,pi/2,om,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus,pi/2);
                subplot(5,6,round(k+(focus-0.8)*60))
                pcolor(cos(phi)*rho,sin(phi)*rho,int1+int2); 
                axis image; shading interp; axis off
                colormap hot
                drawnow
            end
        end
        if kappa<0
            eval(['print -djpeg QDImage-' mint2str(10*abs(kappa),2)]);
        else
            eval(['print -djpeg QDImage' mint2str(10*abs(kappa),2)]);            
        end
    end
end

if 0 % elliptic 2D-emitter
    cnt = 1;
    focus = 1.03;
    for k=1:7
        om = pi/12*(k-1);
        for j=1:11
            kappa = -1+(j-1)/10*2;        
            [int1, rho, phi] = QDSEP(kappa,0,om,0,1.5,0,NA,1.51,1,1,[],0,[],lambda,mag,focus);
            int2 = QDSEP(kappa,pi/2,om,0,1.5,0,NA,1.51,1,1,[],0,[],lambda,mag,focus,pi/2);
            subplot(7,11,cnt)
            pcolor(cos(phi)*rho,sin(phi)*rho,int1+int2); 
            if k==1
                title(mnum2str(kappa,1,2),'fontname','times','fontsize',10);
            end
            if j==1
                h = text(-250,0,int2str(om/pi*180),'fontname','times','fontsize',10);
                set(h,'horizontalalignment','right');
            end
            axis image; shading interp; axis off
            colormap hot
            drawnow
            cnt = cnt + 1;
        end
    end
end

if 0 % elliptic 2D-emitter
    cnt = 1;
    focus = 1.1;
    for k=1:6
        kappa = -1+(k-1)/5;        
        om0 = 0;
        for j=1:11
            if j<=6
                om = pi/10*(j-1);        
            else
                om0 = pi/10*(j-6);        
            end
            [int1, rho, phi] = QDSEP(kappa,0,om,om0,1.5,0,NA,1.51,1,1,[],0,[],lambda,mag,focus);
            int2 = QDSEP(kappa,pi/2,om,om0,1.5,0,NA,1.51,1,1,[],0,[],lambda,mag,focus,pi/2);
            subplot(6,11,cnt)
            pcolor(cos(phi)*rho,sin(phi)*rho,int1+int2); 
            if k==1
                title({['\Omega = ' int2str(om/pi*180) '°'];['\omega_0 = ' int2str(om0/pi*180) '°']},...
                    'fontname','times','fontsize',10,'HorizontalAlignment','center');
            end
            if j==1
                h = text(-250,0,['\kappa = ' mnum2str(abs(kappa),1,2)],'fontname','times','fontsize',10);
                set(h,'horizontalalignment','right');
            end
            axis image; shading interp; axis off
            colormap hot
            drawnow
            cnt = cnt + 1;
        end
    end
end

if 0 % elliptic 2D-emitter
    clear mask
    nn = 20;
    n2 = 2*nn;
    n3 = (2*nn+1)^2;
    [x,y] = meshgrid(-nn:nn,-nn:nn);
    p = angle(x+i*y);
    r = sqrt(x.^2+y.^2);

    cnt = 1;
    focus = 1.15;
    for k=1:7
        om = pi/12*(k-1);
        for j=1:11
            kappa = -1+(j-1)/10*2;        
            [tmp, rho, phi, intc1, ints1] = QDSEP(kappa,0,om,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus);
            [tmp, rho, phi, intc2, ints2] = QDSEP(kappa,pi/2,om,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus,pi/2);
            rho = rho/6.45;
            subplot(7,11,cnt)
            int = interp1(rho,real(intc1(1,:)+intc2(1,:)),r,'cubic');
            for jj=1:4 
                int = int + cos(jj*p).*interp1(rho,real(intc1(jj+1,:)),r,'cubic') + ...
                    sin(jj*p).*interp1(rho,real(ints1(jj,:)),r,'cubic'); 
                int = int + cos(jj*(p+pi/2)).*interp1(rho,real(intc2(jj+1,:)),r,'cubic') + ...
                    sin(jj*(p+pi/2)).*interp1(rho,real(ints2(jj,:)),r,'cubic');                 
            end
            mask(:,:,cnt) = int;
            mim(mask(:,:,cnt));
            if k==1
                title(mnum2str(kappa,1,2),'fontname','times','fontsize',10);
            end
            if j==1
                h = text(-10,21,int2str(om/pi*180),'fontname','times','fontsize',10);
                set(h,'horizontalalignment','right');
            end
            drawnow
            cnt = cnt + 1;
        end
    end
end

if 0 % elliptic 2D-emitter
    clear mask
    nn = 20;
    n2 = 2*nn;
    n3 = (2*nn+1)^2;
    [x,y] = meshgrid(-nn:nn,-nn:nn);
    p = angle(x+i*y);
    r = sqrt(x.^2+y.^2);

    cnt = 1;
    om = 0;
    kappa = 0;
    
    for k=1:7
        focus = 1 + (k-4)*0.01;
        for j=1:11
            NA = 1.25 + (j-1)/10*0.05;        
            [tmp, rho, phi, intc1, ints1] = QDSEP(kappa,0,om,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus);
            [tmp, rho, phi, intc2, ints2] = QDSEP(kappa,pi/2,om,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus,pi/2);
            rho = rho/6.45;
            subplot(7,11,cnt)
            int = interp1(rho,real(intc1(1,:)+intc2(1,:)),r,'cubic');
            for jj=1:4 
                int = int + cos(jj*p).*interp1(rho,real(intc1(jj+1,:)),r,'cubic') + ...
                    sin(jj*p).*interp1(rho,real(ints1(jj,:)),r,'cubic'); 
                int = int + cos(jj*(p+pi/2)).*interp1(rho,real(intc2(jj+1,:)),r,'cubic') + ...
                    sin(jj*(p+pi/2)).*interp1(rho,real(ints2(jj,:)),r,'cubic');                 
            end
            mask(:,:,cnt) = int;
            mim(mask(:,:,cnt));
            if k==1
                mnum2str(focus,1,3)
                title(mnum2str(NA,1,3),'fontname','times','fontsize',10);
            end
            if j==1
                h = text(-10,21,mnum2str(focus,1,3),'fontname','times','fontsize',10);
                set(h,'horizontalalignment','right');
            end
            drawnow
            cnt = cnt + 1;
        end
    end
end


if 0 % elliptic 2D-emitter
    clear mask
    nn = 20;
    n2 = 2*nn;
    n3 = (2*nn+1)^2;
    [x,y] = meshgrid(-nn:nn,-nn:nn);
    p = angle(x+i*y);
    r = sqrt(x.^2+y.^2);

    cnt = 1;
    focus = 1.05;
    kappa = -0.6;
    for k=1:7
        om = pi/12*(k-1);
        for j=1:11
            om0 = (j-1)/10*pi/2;        
            [tmp, rho, phi, intc1, ints1] = QDSEP(kappa,0,om,om0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus);
            [tmp, rho, phi, intc2, ints2] = QDSEP(kappa,pi/2,om,om0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus,pi/2);
            rho = rho/6.45;
            subplot(7,11,cnt)
            int = interp1(rho,real(intc1(1,:)+intc2(1,:)),r,'cubic');
            for jj=1:4 
                int = int + cos(jj*p).*interp1(rho,real(intc1(jj+1,:)),r,'cubic') + ...
                    sin(jj*p).*interp1(rho,real(ints1(jj,:)),r,'cubic'); 
                int = int + cos(jj*(p+pi/2)).*interp1(rho,real(intc2(jj+1,:)),r,'cubic') + ...
                    sin(jj*(p+pi/2)).*interp1(rho,real(ints2(jj,:)),r,'cubic');                 
            end
            mask(:,:,cnt) = int;
            mim(mask(:,:,cnt));
            if k==1
                title(int2str(om0/pi*180),'fontname','times','fontsize',10);
            end
            if j==1
                h = text(-10,21,int2str(om/pi*180),'fontname','times','fontsize',10);
                set(h,'horizontalalignment','right');
            end
            drawnow
            cnt = cnt + 1;
        end
    end
end

if 0 % 3D-isotropic emitter
    for k=1:length(focus)
        [intx inty int, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
            SEPDipole([0 5], 0, NA, 1.51, 1., 1., [], 0, [], lambda, mag, focus(k));
        col = ones(size(phi));
        int = real(col*(fxx0.*conj(byx0)+fxz.*conj(byz)+fxx2.*conj(byx2)));
        subplot(5,8,k)
        pcolor(cos(phi)*rho,sin(phi)*rho,int);axis image; shading interp; axis off
        title(num2str(focus(k)),'fontsize',12);
        colormap(gray)
        drawnow
    end
end

if 0 % 3D-isotropic emitter with excitation weighting
    for k=1:length(focus)
        [intx inty int, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
            SEPDipole([0 5], 0, NA, 1.51, 1., 1., [], 0, [], lambda, mag, focus(k));
        col = ones(size(phi));
        int = real(col*(6*fxx0.*conj(byx0)+2*fxz.*conj(byz)+6*fxx2.*conj(byx2))-3*cos(2*phi)*(fxx2.*conj(byx0)+fxx0.*conj(byx2)));
        subplot(5,8,k)
        pcolor(cos(phi)*rho,sin(phi)*rho,int);axis image; shading interp; axis off
        title(num2str(focus(k)),'fontsize',12);
        colormap(gray)
        drawnow
    end
end

if 0 % quadrupole
    for k=1:length(focus)
        [intx inty int, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
            SEPQuadrupole([0 3], 0, NA, 1.51, 1., 1., [], 0, [], lambda, mag, focus(k), [], 1);
        subplot(5,8,k)
        pcolor(cos(phi)*rho,sin(phi)*rho,intx+inty);axis image; shading interp; axis off
        title(num2str(focus(k)),'fontsize',12);
        colormap(gray)
        drawnow
    end
end
