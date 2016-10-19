cd 'c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Pi Ethylenglycol\Lifetimes\'

clear all
close all

nglass = 1.52; % glass
nglycol = 1.43;
dagbottom = 0.036;
dagtop = 0.072;

load SilverDmitri
load metals
lamem = 0.535;
lamex = 0.47;

nagex = interp1(SilverDmitri(:,1),SilverDmitri(:,2),1e3*lamex);
nagem = interp1(SilverDmitri(:,1),SilverDmitri(:,2),1e3*lamem);
nauex = interp1(wavelength,gold,1e3*lamex);
nauem = interp1(wavelength,gold,1e3*lamem);

rhofield = [0 1];
NA = 1.25;
fd = 164.5e3/100;
over = 3e3;
focpos = 0;
mag = 100;
av = 75;
zpin = 0;
atf = [];
kappa = 1;
d = 0.12;

n0 = [nglass nagex nauex]; d0 = [dagbottom 0.004];
n = nglycol;
n1 = [nauex nagex nglass]; d1 = [0.004 dagtop];

resolution = [20 lamex/(d/100)];
exc = GaussExc(rhofield, [0, d], NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);

n0 = [nglass nagem nauem];
n1 = [nauem nagem nglass];

mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, 1);
tau = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, 2);
tau1 = 4/3*n*squeeze(sum(sum(exc.rho.*(tau.volx(:,:,1)+tau.voly(:,:,1)))))/squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))));
brightness1 = squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))))*diff(exc.z(1,1:2));

mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, 0, 1);
lifetime = 4*n./(2*mdf.qp+mdf.qv)';
weight = squeeze(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1))))';
tau0 = weight'*lifetime/sum(weight);
brightness0 = squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))))*diff(exc.z(1,1:2));
