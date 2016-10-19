fnames = {'c:\Users\Public\Transfer\DSOM_beads\5to25.mat',...
    'c:\Users\Public\Transfer\DSOM_beads\5to50.mat',...
    'c:\Users\Public\Transfer\DSOM_beads\5to100.mat'};

for k=1:length(fnames)
    load(fnames{k});
    j=1;
    [mx(k,j), my(k,j), wx(k,j), wy(k,j), amp(k,j)] = Gauss2D([], [], low, [], 1);
    j=2;
    [mx(k,j), my(k,j), wx(k,j), wy(k,j), amp(k,j)] = Gauss2D([], [], high, [], 1);
    j=3;
    [mx(k,j), my(k,j), wx(k,j), wy(k,j), amp(k,j)] = Gauss2D([], [], test, [], 1);
end