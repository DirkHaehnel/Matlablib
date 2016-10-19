resolution = [60 20];
lamex = 0.635;

rhofield = [0 0.8];
zfield = [-1 4];
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
modres = FCS(exc.rho,exc.z,mdf.volx,mdf.voly);

dz = exc.z(1,2)-exc.z(1,1);
save CuttedFCSRes exc mdf modres
for j=1:size(exc.z,2)-1
    tic
    [auto3d(:,j), autotime] = FCS(exc.rho(:,1+j:end),exc.z(:,1+j:end)-exc.z(1,1+j-1)-dz/2,mdf.volx(:,1+j:end,:),mdf.voly(:,1+j:end,:));
    auto2d(:,j) = FCS(exc.rho(:,j),[],mdf.volx(:,j,:),mdf.voly(:,j,:));
    save CuttedFCSRes auto3d auto2d autotime -append
end

for j=1:size(exc.z,2)-1
    if sum(auto3d(:,j))>0 && all(isfinite(auto3d(:,j)))
        tau(j) = interp1(auto3d(diff(auto3d(:,j))<0,j),autotime(diff(auto3d(:,j))<0),0.5*(1+auto3d(end,j)),'cubic');
    else
        tau(j) = 0;
    end
end

tau0 = interp1(modres(diff(modres)<0),autotime(diff(modres)<0),0.5*(1+modres(end)),'cubic');
save CuttedFCSRes tau tau0 -append

