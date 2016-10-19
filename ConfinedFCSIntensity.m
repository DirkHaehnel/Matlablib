radiusv = [0.1 0.2 0.5 1 2 3];
rv = (0:60)/60;
pv = 0:pi/36:pi/2;

resolution = 60;
lamex = 0.635;
 
rhofield = [-lamex/resolution/2 0.5+lamex/resolution];
zfield = [-lamex/resolution/2 1+lamex/resolution];
NA = 1.2;
fd = 5e3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = [];
d1 = [];
over = 5e3;
focpos = 0;
atf = [];
ring = '';
lamem = 0.67;
mag = 60;
av = 10; %pinhole radius!!!
zpin = 0;

exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf);

phiv = (0.5:360)/180*pi;
[zz,rr,pp] = meshgrid([-fliplr(exc.z(1,2:end)) exc.z(1,:)],exc.rho(:,1)',phiv);
rc = rr.*cos(pp); rs2 = (rr.*sin(pp)).^2;
dist = rr.*repmat([fliplr(mdf.volx(:,2:end,1)+mdf.voly(:,2:end,1)) mdf.volx(:,:,1)+mdf.voly(:,:,1)],[1 1 length(phiv)]);

for jr=1:length(radiusv)
    radius = radiusv(jr); % cell radius in mum
    for kr=1:length(rv)
        for jp=1:length(pv)
            x0 = (radius+2)*rv(kr)*sin(pv(jp));
            z0 = (radius+2)*rv(kr)*cos(pv(jp));
            ind = (x0+rc).^2 + rs2 + (z0+zz).^2<=radius^2;
            if sum(ind(:))>0
                int(kr,jp,jr) = sum(dist(ind));
            else
                int(kr,jp,jr) = 0;
            end
        end
        kr
    end
end

int0 = sum(dist(:));

save ConfinedFCSInt exc mdf radiusv pv rv int

for j=1:length(radiusv)
    tmp = int(:,:,j)/int0;
    pcolor((2+radiusv(j))*rv'*sin(pv),(2+radiusv(j))*rv'*cos(pv),tmp);
    colormap hot
    hold on;
    pcolor(-(2+radiusv(j))*rv'*sin(pv),(2+radiusv(j))*rv'*cos(pv),tmp);
    pcolor((2+radiusv(j))*rv'*sin(pv),-(2+radiusv(j))*rv'*cos(pv),tmp);
    pcolor(-(2+radiusv(j))*rv'*sin(pv),-(2+radiusv(j))*rv'*cos(pv),tmp);
    axis image;
    shading interp
    plot(radiusv(j)*sin(0:pi/100:2*pi),radiusv(j)*cos(0:pi/100:2*pi),'y'); hold off
    colorbar
    patch([-0.5 0.5 0.5 -0.5],[-radiusv(j)-1.55 -radiusv(j)-1.55 -radiusv(j)-1.45 -radiusv(j)-1.45],'w')
    axis off
    eval(['print -dpng -r300 ConfinedInt' mint2str(10*radiusv(j),2)]);
end



