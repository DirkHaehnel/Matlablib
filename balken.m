function balken(x,s)

r = 0.12*(x(end)-x(1));
x = x(:)';

if nargin>1 & s=='v'
    pcolor([0;r]*ones(1,length(x)),[x;x],[x;x]);
    set(gca,'xtick',[],'yaxislocation','right','fontsize',20)
else
    pcolor([x;x],[0;r]*ones(1,length(x)),[x;x]);
    set(gca,'ytick',[],'fontsize',28)
end

axis image
colormap hot



