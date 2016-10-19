% Emission rates for the radiation of a
% dipole interacting with a cylindrical structure

close all
clear

load metals
lamem = 570;
nh2o = 1.33;
nal = silver(wavelength==lamem);
nglass = 1.52;
inner_radius = 75;
metal_thickness = 5e2;
%dipole_pos_vec = [1e-2 5:2.5:70 75-1e-2];
dipole_pos_vec = [0:5:70 71:74];
dipole_pos_vec(1) = 1e-2;

n0 = nh2o;
n = nh2o;
n1 = [nal nglass];
drho = 0;
%nmax = 35;
nmax = 50;

for jr = 1:length(dipole_pos_vec)
    dipole_pos = dipole_pos_vec(jr);
    d0 = dipole_pos/lamem*2*pi;
    d = (inner_radius - dipole_pos)/lamem*2*pi;
    d1 = metal_thickness/lamem*2*pi;

    [lr(jr),lf(jr),lz(jr),qrout(jr),qfout(jr),qzout(jr),qrin(jr),qfin(jr),qzin(jr)] = LifetimeC1(nmax,drho,n0,n,n1,d0,d,d1);

end

%save LifetimeSilverCylinder75-10Ag

ri = 0:0.1:75; 
semilogy(ri,exp(interp1(dipole_pos_vec,log(qrout/(4/3*n)),ri,'cubic')),ri,exp(interp1(dipole_pos_vec,log(qfout/(4/3*n)),ri,'cubic')),ri,exp(interp1(dipole_pos_vec,log(qzout/(4/3*n)),ri,'cubic')))
ax = axis;
axis([0 max(dipole_pos_vec) ax(3:4)]);
xlabel('dipole position (nm)');
ylabel('rel. emission rate');
legend({'\rho-dipole','\phi-dipole','\itz\rm-dipole'},2)

figure
ri = 0:0.1:75; 
semilogy(ri,exp(interp1(dipole_pos_vec,log(lr./qrout),ri,'cubic')),ri,exp(interp1(dipole_pos_vec,log(lf./qfout),ri,'cubic')),ri,exp(interp1(dipole_pos_vec,log(lz./qzout),ri,'cubic')))
ax = axis;
axis([0 max(dipole_pos_vec) ax(3:4)]);
xlabel('dipole position (nm)');
ylabel('transmitted fluorescence');
legend({'\rho-dipole','\phi-dipole','\itz\rm-dipole'})

