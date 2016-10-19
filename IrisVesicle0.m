% IrisVesicle

path = 'D:\Doc\Iris\Vesicles\';
name = 'wtmRecX647_253C.mat';

global pd
pd = [0.01 1e5];
FCSFit([path name]);
pd0 = pd(1);

pd = [1e-4; pd];

files = dir([path '*mgml*.mat'])
for j=1:length(files)
    pos = findstr(files(j).name,'mgml');
    t(j) = str2num(files(j).name(1:pos-1));
    [dc v conc w0 a0 triplet c err] = FCSFit([path files(j).name],[],1,[],[],[[0 pd0 0];inf pd0 inf]);
    dc1(j) = dc(1); dc2(j) = dc(2);
    c1(:,j) = c(:,1); c2(:,j) = c(:,2);
end
    
return

bar(t([7 1:2:5]),[dc1([7 1:2:5]);dc1([8 2:2:6])]');legend({'Ca','EGTA'})
xlabel('mgml')
ylabel('slow diff. coef.')

figure
bar(t([7 1:2:5]),[mean([c1(2,[7 1:2:5])./sum(c1(2:3,[7 1:2:5]));c2(2,[7 1:2:5])./sum(c2(2:3,[7 1:2:5]))]);mean([c1(2,[8 2:2:6])./sum(c1(2:3,[8 2:2:6]));c2(2,[8 2:2:6])./sum(c2(2:3,[8 2:2:6]))])]');legend({'Ca','EGTA'})
xlabel('mgml')
ylabel('\chi_{slow}')
