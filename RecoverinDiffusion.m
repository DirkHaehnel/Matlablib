% Recoverin

tcpath = 'D:\Joerg\Doc\Iris\'
global pd

clear dc v conc w0 a0 triplet
pd=[0.01 10 10 10];
j = 1;
[dc{j} v{j} conc{j} w0{j} a0{j} triplet{j}] = FCSFit('D:\Joerg\Doc\Iris\RecX647_Lösung_74um_12uW_26,3C_neueJustage.mat',[400 150],3);

[y,t] = FCSCrossRead('D:\Joerg\Doc\Fcs\2Focus\GUV\2006-08-29\LoesungRecX647_85um.mat');

pd = [1e-8/dc{1}; 1e-8/dc{1}*10; triplet{1}];

p = simplex('GaussFCS',[w0{1} a0{1} 0 w0{1} w0{1}],[w0{1} a0{1} -inf 0 0],[w0{1} a0{1} inf inf inf],[],[],100e3/58,[640 670]/1.33,...
    403,1e6*t,sum(y(:,1:2,:),3),sum(y(:,3,:),3),3,[[1e-8/dc{1} 1e-8/dc{1} 0 0 0];[1e-8/dc{1} inf 100 100 100]]);