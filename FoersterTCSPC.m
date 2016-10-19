% fitting data from Ute Foerster

load c:\Joerg\Doc\Nanophotonics&Optics\FoersterWachtveitl\Foersterdata.mat

ind = irf(:,1)>0.01 & irf(:,1)<0.06;

[c, offset, A, tau, dc, dtau, irs, zz, t, chi] = fluofit(irf(ind,2), y(ind,2), 1e3, 1);

