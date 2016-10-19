% DefocHeight

focposv = -(-1:0.05:-0.05);
dv = (0:50:250)/1e3;

rhov = [0 2];
z = 0;
NA = 1.242; 
n1 = 1.52;
n = 1.49;
n2 = 1.0;
d1 = 1;
d2 = [];
lambda = 0.57;
mag = 100;

for jd = 1:length(dv)
    for jf=1:length(focposv)
        [intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
            SEPDipole(rhov, z, NA, [n1 n], n, n2, d1, dv(jd), d2, lambda, mag, focposv(jf));
        subplot(4,5,jf);
        mpcolor(cos(phi)*rho,sin(phi)*rho,intx+inty);
        axis off
        title(num2str(focposv(jf)),'fontsize',12);
    end
    suptitle([num2str(dv(jd)) ' nm PMMA layer']);
    eval(['print -dpng -r300 PMMA_' mint2str(1e3*dv(jd),3) 'nm']);
end

