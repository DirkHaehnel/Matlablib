clear
%clc

% Defining range of investigated refractive index.
n1_anfang = 1.333;
n1_ende = 1.400;

% Simulation steps, 2 at least.
num_of_runs_ref = 15;

NA = 1.2;
mag = 60;
tubelens = 164.5e3;
fd = tubelens/mag;
d0 = [];
d = [];
d1 = [];
lamex = 0.5125;
over = 5e3;
focpos = 200;
foc_shift = 7;
xy_size = 1;
lamem = 0.570;
av = 35; %pinhole radius!!!
zpin = 0;
tsh = 1/exp(2);
rhofield = [0 1];
zfield = [focpos-foc_shift focpos+foc_shift];
n0 = 1.333;

for j = 1:num_of_runs_ref

    n1 = n1_anfang*exp(log(n1_ende/n1_anfang)*(j-1)/(num_of_runs_ref-1));
    n = n1;

    % determine fcous position
    shft = GaussExc(0, [150 250], NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos);
    [ind,ind] = max(abs(shft.fxc(1,:,1)));
    shft = shft.z(1,ind)-focpos;
    
    % calculate excitation and detection
    exc = GaussExc(rhofield, shft+zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos);
    mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin);

    % draw image
    FocusImage3D(exc.rho,exc.z,mdf.volx+mdf.voly)
    text(-3,3,shft+focpos,['\itn\rm_1 = ' mnum2str(n1,1,3)]);
    axis([-xy_size xy_size -xy_size xy_size shft+focpos-foc_shift shft+focpos+foc_shift])
    view([70 20])
    axis vis3d
    drawnow;
    eval(['print -dpng -r300 RefIndexMDF' mint2str(j,2)]);
    if j<num_of_runs_ref
        eval(['print -dpng -r300 RefIndexMDF' mint2str(2*num_of_runs_ref-j,2)]);
    end

end
