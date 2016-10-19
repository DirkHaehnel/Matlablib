% RecoverinMembraneAnalysis

% Fitting of solution 2fFCS for getting the diffusion of free Recoverin
global pd
pd = [0.02 10 10 10];
%[dc v conc w0 a0 triplet c err] = FCSFit('D:\Joerg\Doc\Fcs\2Focus\GUV\2006-08-29\LoesungRecX647_85um.mat');
[y,t] = FCSCrossRead('D:\Joerg\Doc\Fcs\2Focus\GUV\2006-08-07 RecX647 DOPC-GUVs\RecX647_DOPC_P11_2uW_75,0um_25,6C.mat');
data.y=y(t>5e-6,:); data.t=t(t>5e-6);
[dc v conc w0 a0 triplet c err] = FCSFit(data);

% Fixing the diffusion of free Recoverin
dc = 0.98e-6;

% Fitting membrane + solution
% meaning of p: 
% p(1) - beam waist w0 for solution fitting
% p(2) - pinhole function parameter a0 for solution fitting
% p(3) - membrane position with respect to beam waist
% p(4) - beam waist within membrane for membrane fitting

[y,t] = FCSCrossRead('D:\Joerg\Doc\Fcs\2Focus\GUV\2006-08-07 RecX647 DOPC-GUVs\RecX647_DOPC_P11_2uW_65,5um_25,6C.mat');

% Display of all curves
semilogx(t,squeeze(y(:,1,:))./(ones(size(y,1),1)*squeeze(y(end,1,:))'))

% separate plot of 3rd curve: semilogx(t,squeeze(y(:,:,3)))
% If you would like to delete the 3rd and 5th curves: y(:,:,[3 5]) = [];

% Start of fitting
pd = [1e-8/dc 0.1 40];
p = [350 150 -100 350];
for j=1:21
    p(3) = -1100 + j*100;
    close; p = Simplex('GaussFcs',p,[],[],[],[],100e3/58,[640 670]/1.33,400,1e6*t,sum(y(:,1:2,:),3),sum(y(:,3,:),3),1,[[1e-8./dc 1e-8./dc 0];[1e-8./dc inf inf]])
    res.p(:,j) = p;
    res.pd(:,j) = pd;
    eval(['print -dpng -r300 Fit' mint2str(j,2)])
end
% results
% p = [534   80 -174  395]; 
% pd = [0.0105 0.0863 46]

[y,t] = FCSCrossRead('D:\Joerg\Doc\Fcs\2Focus\GUV\2006-08-29\BilayerRecX647_70um.mat');
pd = [1e-8/dc 0.1 40];
close; p=Simplex('GaussFcs',p,[],[],[],[],100e3/58,[640 670]/1.33,373,1e6*t,sum(y(:,1:2,:),3),sum(y(:,3,:),3),1,[[1e-8./dc 0 0];[1e-8./dc inf inf]])
% results
% p = [537   82 -184  405]; 
% pd = [0.0105 0.0947 53]

[y,t] = FCSCrossRead('D:\Joerg\Doc\Fcs\2Focus\GUV\2006-08-29\BilayerRecX647_70,5um.mat');
pd = [1e-8/dc 0.1 40];
close; p=Simplex('GaussFcs',p,[],[],[],[],100e3/58,[640 670]/1.33,373,1e6*t,sum(y(:,1:2,:),3),sum(y(:,3,:),3),1,[[1e-8./dc 0 0];[1e-8./dc inf inf]])
% results
% p = [549  80 -184  476]; 
% pd = [0.0105 0.1126 36]

[y,t] = FCSCrossRead('D:\Joerg\Doc\Fcs\2Focus\GUV\2006-08-29\BilayerRecX647_71um.mat');
pd = [1e-8/dc 0.1 40];
close; p=Simplex('GaussFcs',p,[],[],[],[],100e3/58,[640 670]/1.33,373,1e6*t,sum(y(:,1:2,:),3),sum(y(:,3,:),3),1,[[1e-8./dc 0 0];[1e-8./dc inf inf]])
% results
% p = [581 80 -184 471]; 
% pd = [0.0105 0.1622 35]




