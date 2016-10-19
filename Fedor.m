% program Fedor Jelezko

rhofield = [0 1.5];
zfield = [0 15];
NA = 1.65;
n = 2.418; % diamond
ng = 1.78;
mag = 100;
over = [NA*200/mag*1e3 0 inf];
zpin = 0; 
atf = []; %1.48;
kappa = 1;
lt = 1;
sat = 0;
lamex = 0.514;
lamem = 0.57;
av = 50/2;
resolution = 50;

fv = 0:0.2:10;
for k=1:length(fv)
    focpos = fv(k);
    [volx, voly, rho, z, h0, h1, h2, qv, qp] = ...
        MDF(rhofield, zfield, NA, ng, n, n, [], 0, [], lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution);
    feld = [volx(:,:,1)+voly(:,:,1)];
    mpc(rho,z,feld,'lr');
    xlabel('\rho [\mum]');
    ylabel('\itz\rm [\mum]');
    title(['foc. pos. ' mnum2str(focpos,2,1) ' \mum']);
    eval(['print -dpng -r110 FedorMDF' mint2str(k,2)])
    mpc(rho,z,abs(h0).^2+abs(h2).^2+abs(h1).^2/2,'lr');
    xlabel('\rho [\mum]');
    ylabel('\itz\rm [\mum]');
    title(['foc. pos. ' mnum2str(focpos,2,1) ' \mum']);
    eval(['print -dpng -r110 FedorExc' mint2str(k,2)])
end


