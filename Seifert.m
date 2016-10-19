% calculation of mplecule-to-surface distance sensitivity of defoc imaging

rhov = [0 2]; 
zv = 0:0.01:0.2;
NA = 1.4;
n1 = 1.51;
n = 1.333;
n2 = 1.333;
d1 = [];
d2 = [];
lambda = 0.514;
mag = 100;
focpos = 1;
atf = [];

for j=1:length(zv)
    [intx, inty, intz, rho, phi] = SEPDipole(rhov, zv(j), NA, n1, n, n2, d1, zv(j), d2, lambda, mag, focpos);
    subplot(3,7,j)
    %pcolor(cos(phi)*rho,sin(phi)*rho,intx+inty); 
    pcolor(cos(phi)*rho,sin(phi)*rho,intz); 
    colormap hot; 
    shading interp; 
    axis image;  
    axis off
    title([int2str(zv(j)*1e3) ' nm'])
    drawnow
end
    
