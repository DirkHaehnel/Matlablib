% program AxiconExc

rhofield = [0 2];
zfield = [-3 3];
NA = 1.2;
fd = 3e3;
n0 = 1.33;
n = 1.33;
n1 = 1.33;
d0 = [];
d = 0;
d1 = [];
lamex = 0.514;
over = 5e3;
focpos = 0;
atf = [];
resolution = [25 5];

tsh = 10:30;

for j=1:length(tsh)
    ring = ['exp(i*2*pi*' num2str(tsh(j)) '*rad)']; 
    zfield = [-10 10] + tsh(j)*0.5;
    [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring);
    FocusImage2D(rho,z,cat(4,cat(3,fxc,fxs),cat(3,fyc,fys),cat(3,fzc,fzs)))
    ax = axis;
    text(0,ax(3)+0.9*diff(ax(3:4)),['\Delta\psi_{max} = ' mnum2str(tsh(j),1,2) '\cdot2\cdot\pi'],'Color','w')
    eval(['print -dpng tap' mint2str(j,2)]);
end

