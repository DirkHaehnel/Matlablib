% example for flux measurement fitting

name = 'c:\Joerg\Doc\Anastasia\Channel\beads 25C pressure 10mbar.mat';

global pd
pd = 10;
p = [350 150 -10 10];
[dc v conc w0 a0 triplet c velo] = FCSFit(name, p);
% refinement of fit:
for j=1:3
    [dc v conc w0 a0 triplet c velo err] = FCSFit(name,[w0 a0 velo(1:2)']);
end


