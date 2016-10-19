function res = CurveSelect

if nargin==0
    h = flipud(get(gca,'children'));
end
for j=1:length(h) 
    x{j} = get(h(j),'xdata');
    y{j} = get(h(j),'ydata');
end

[a,b] = ginput;

for j=1:length(a)
    tst = zeros(1,length(x));
    for k=1:length(x)
        tst(k) = abs(b(j)-interp1(x{k},y{k},a(j)));
    end
    [tmp,res(j)] = min(tst);
end

