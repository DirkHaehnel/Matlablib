% program for calculating detection volumina and overlap for FRET measurements

lambda = [488 511.3 530 589.3 595.6 623.4 635 671 770 820 880 902.2];
water = [1.33611 1.33525 1.33447 1.33236 1.33217 1.33138 1.33108 1.33022 1.32828 1.32741 1.32641 1.32606];
urea = [1.40418 1.40293 1.40181 1.39875 1.39849 1.39741 1.39698 1.39583 1.39318 1.39208 1.39098 1.39063];
gdcl = [1.42491 1.42345 1.42213 1.41858 1.41827 1.41703 1.41654 1.41522 1.41225 1.41103 1.4098 1.40941];
crys = [1.3838 1.38268 1.38167 1.37892 1.37868 1.3777 1.37732 1.37625 1.37386 1.37287 1.3718 1.37141]; % gamma crystallin

solution = {'water','urea','gdcl','crys'};

over = 2e3;
focpos = 100;
rhofield = [0 10];
zfield = focpos + [-15 15];
NA = 1.2;
d0 = [];
d = 0;
d1 = [];
lamex = 0.488;
lamem = [0.53 0.67];
n0ex = interp1(lambda,water,lamex*1e3);
n0em = interp1(lambda,water,lamem*1e3);
mag = 60;
av = 50;
zpin = 0;
atf = [];
kappa = 1;
lt = 0;
fd = 3e3;
resolution = [30 10];

for sol=1:length(solution)
    eval(['nex = interp1(lambda,' solution{sol} ',lamex*1e3);']);
    eval(['nem = interp1(lambda,' solution{sol} ',lamem*1e3);']);

    % determine fcous position
    shft = GaussExc(0, focpos + [-50 50], NA, fd, n0ex, nex, nex, d0, d, d1, lamex, over, focpos);
    [ind,ind] = max(abs(shft.fxc(1,:,1)));
    shft = shft.z(1,ind)-focpos;
    
    exc = GaussExc(rhofield, shft + zfield, NA, fd, n0ex, nex, nex, d0, d, d1, lamex, over, focpos, atf, resolution);
    mdf(1) = GaussExc2MDF(exc, NA, n0em(1), nem(1), nem(1), focpos, lamem(1), mag, av, zpin, atf, kappa, lt);
    mdf(2) = GaussExc2MDF(exc, NA, n0em(2), nem(2), nem(2), focpos, lamem(2), mag, av, zpin, atf, kappa, lt);

    rho = exc.rho(:,1);
    z = exc.z;
    maxm = size(exc.fxs,3);

    subplot(121)
    FocusImage2D(exc.rho,exc.z,mdf(1).volx+mdf(1).voly)
    axis([-2 2 shft+focpos+[-4 4]])
    text(0,shft+focpos+3.5,['\lambda_{\item\rm} = ' mint2str(1e3*lamem(1),3)],'color','y','horizontalalignment','center')
    subplot(122)
    FocusImage2D(exc.rho,exc.z,mdf(2).volx+mdf(2).voly)
    axis([-2 2 shft+focpos+[-4 4]])
    text(0,shft+focpos+3.5,['\lambda_{\item\rm} = ' mint2str(1e3*lamem(2),3)],'color','y','horizontalalignment','center')
    eval(['print -dpng -r300 Khambiz_' solution{sol}])
    
    v1(sol)  = DetectionVolume(exc.rho,exc.z,mdf(1).volx+mdf(1).voly);
    v2(sol) = DetectionVolume(exc.rho,exc.z,mdf(2).volx+mdf(2).voly);
    vc(sol) = DetectionVolume(exc.rho,exc.z,mdf(1).volx+mdf(1).voly,mdf(2).volx+mdf(2).voly);
end

