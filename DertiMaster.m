radv = 2.6:0.2:4;
disv = 0.19:0.01:0.21
% load D:\Doc\Dertinger\Methanol\atto655_d-methanol_4e6_10uW_2.mat
load D:\MATLAB\Cy5_e-10_h2o_25C_005uW.mat 

fcs1=[res.auto(:,1,2)/res.auto(end,1,2) res.auto(:,2,1)/res.auto(end,2,1)]-1;
fcs2=[res.auto(:,3,4)/res.auto(end,3,4) res.auto(:,4,3)/res.auto(end,4,3)]-1;
xfcs1=(res.auto(:,1,4)/res.auto(end,1,4)+res.auto(:,3,2)/res.auto(end,3,2))-2;
xfcs2=(res.auto(:,2,3)/res.auto(end,2,3)+res.auto(:,4,1)/res.auto(end,4,1))-2;
xfcs = [xfcs1 xfcs2];


for jdis=1:length(disv)
    focpos = [15 disv(jdis) 0 0 0];
    for jrad=1:length(radv)
        over = [0 radv(jrad)*1e3/2];
        TwoFocusFCS

        p1=simplex('AffineExpFit',[1e-6 1e-5 1e-5 1e-5],zeros(1,4),[],[],[],res.autotime,[fcs1(:,1) fcs2(:,1)],autotime,modres, xfcs1);
        for k=1:2
            p1=simplex('AffineExpFit',p1,0*p1,[],[],[],res.autotime,[fcs1(:,1) fcs2(:,1)],autotime,modres, xfcs1);
        end
        
        p2=simplex('AffineExpFit',[1e-6 1e-5 1e-5 1e-5],zeros(1,4),[],[],[],res.autotime,[fcs1(:,2) fcs2(:,2)],autotime,modres, xfcs2);
        for k=1:2
            p2=simplex('AffineExpFit',p2,0*p2,[],[],[],res.autotime,[fcs1(:,2) fcs2(:,2)],autotime,modres, xfcs2);
        end

        err1(jrad,jdis)=AffineExpFit(p1,res.autotime,[fcs1(:,1) fcs2(:,1)],autotime, modres, xfcs1);
        err2(jrad,jdis)=AffineExpFit(p2,res.autotime,[fcs1(:,2) fcs2(:,2)],autotime, modres, xfcs2);

        dres(:,jrad,jdis) = 5e-5*1e-8./[p1(1);p2(1)];
    end
end
