lamex = 0.491;
lamem = 0.525;
n0 = 1.33;
n = 1.33;
n1 = 1.33;
d0 = [];
d = 0.1;
d1 = [];
NA = 0;
mag = 60;
fd = 164.5e3/mag;
av = 10;
resolution = [lamex/0.02 lamex/0.02]; % spatial grid resolution for MDF calculation 
focpos = 0; % position of focus in mum
over = [0 0];
atf = [];
zpin = 0;
kappa = 0;
sat = 0;

% calculating the detection efficiency distribution
exc.NA = NA;
exc.fd = fd;
exc.n0 = n0;
exc.n = n;
exc.n1 = n1;
exc.d0 = d0/2/pi*lamex;
exc.d = d/2/pi*lamex;
exc.d1 = d1/2/pi*lamex;
exc.focpos = focpos;
exc.atf = atf;
exc.ring = [];
exc.maxm = 0;

drho = lamex/resolution(1);
dz = lamex/resolution(2);
rhov = 0;
zv = 0;
[exc.rho, exc.z] = ndgrid(rhov, zv);
exc.fxc = ones(size(exc.rho));
exc.fyc = zeros(size(exc.rho));
exc.fzc = zeros(size(exc.rho));
exc.fxs = [];
exc.fys = [];
exc.fzs = [];

[ux, uy, pp, po, op] = GaussExc2SatRotoFCS(exc, lamem, mag, av, zpin, atf, sat);

