% calculate PSF inside cell

if 1

    over = 5e3;
    cov = 0;
    rhofield = [0 1];
    NA = 1.2;
    n0 = 1.333;
    ncyto = 1.36;
    nnucl = 1.4;
    d0 = [];
    dcyto1 = 4;
    dcyto2 = 4;
    dnucl = 2;
    lamex = 0.514;
    lamem = 0.55;
    mag = 40;
    av = 45;
    zpin = 0e3;
    atf = [];
    lt = [];
    kappa = 1;
    fd = 164.5e3/mag;
    resolution = 2*[20 20];

    focposv = 0:2:10;
    
    for j=1:length(focposv)
        focpos = focposv(j);
        exc1 = GaussExc(rhofield, [0 dcyto1], NA, fd, n0, ncyto, [nnucl ncyto n0], [], dcyto1, [dnucl dcyto2], lamex, over, focpos, atf, resolution);
        exc2 = GaussExc(rhofield, [0 dnucl], NA, fd, [n0 ncyto], nnucl, [ncyto n0], [dcyto1], dnucl, dcyto2, lamex, over, focpos, atf, resolution);
        exc3 = GaussExc(rhofield, [0 dcyto2], NA, fd, [n0, ncyto, nnucl], ncyto, n0, [dcyto1, dnucl], dcyto2, [], lamex, over, focpos, atf, resolution);
       
        mdf1 = GaussExc2MDF(exc1, NA, n0, ncyto, [nnucl ncyto n0], focpos, lamem, mag, av, zpin);
        mdf2 = GaussExc2MDF(exc2, NA, [n0 ncyto], nnucl, [ncyto n0], focpos, lamem, mag, av, zpin);
        mdf3 = GaussExc2MDF(exc3, NA, [n0, ncyto, nnucl], ncyto, n0, focpos, lamem, mag, av, zpin);
       
        subplot(1,length(focposv),j)
        FocusImage2D([exc1.rho exc2.rho exc3.rho],[exc1.z dcyto1+exc2.z dcyto1+dnucl+exc3.z],cat(2,mdf1.volx+mdf1.voly,mdf2.volx+mdf2.voly,mdf3.volx+mdf3.voly));
    end

end