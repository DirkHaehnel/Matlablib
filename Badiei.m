if 1 % Opt. Express 2006
    
    close all
    clear all
    
    rhov = [0 0.8];
    z = 0;
    NA = 1.4;
    n1 = 1.52;
    n = 1.33;
    n2 = 1.33;
    d1 = [];
    d = 0;
    d2 = [];
    lambda = 0.55;
    mag = 100;
    focpos = 0.0;
    atf = [];
    ring = [];
    
    thetav = (0:10:90)/180*pi;
    
    for j=1:length(thetav)
        orient = [thetav(j) 0];
        
        [int, ~, ~, rho, phi] = SEPDipole(rhov, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, orient);
        rho = rho/mag;
        subplot(2,5,j)
        pcolor(cos(phi)*rho,sin(phi)*rho,int); axis image; shading interp
    end
    colormap hot
    
end

if 0 % JOSA B 2003
    
    close all
    
    rhov = [0 4];
    z = 0;
    NA = 1.4;
    n1 = 1.51;
    n = 1.0;
    n2 = 1.0;
    d1 = [];
    d = 0;
    d2 = [];
    lambda = 0.647;
    mag = 60;
    atf = [];
    ring = [];
    
    thetav = fliplr(0:45:90)/180*pi;
    focusv = [0.6, 1.2];
    cnt = 1;
    for jf=1:2
        focpos = focusv(jf);
        for jt=1:3
            orient = [thetav(jt) pi/2];
            [int, ~, ~, rho, phi] = SEPDipole(rhov, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ring, orient);
            subplot(2,3,cnt)
            pcolor(cos(phi)*rho,sin(phi)*rho,int); axis image; shading interp
            cnt = cnt+1;
        end
    end
    colormap hot
    
end

