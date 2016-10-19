function ind = PickData()

tst = get(gca,'children');
j = 1;
while strcmp(get(tst(j),'Marker'),'none')
    j = j+1;
end

xdata = get(tst(j),'XData');
ydata = get(tst(j),'YData');

[x,y] = ginput;

v = 1:length(xdata);
for k=1:length(x)
    tst = (xdata-x(k)).^2+(ydata-y(k)).^2;
    ind(k) = v(tst==min(tst));
end
