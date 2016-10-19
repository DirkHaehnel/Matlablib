function colorize(h)

if nargin==0
    h = get(gca,'children');
end
for j=1:length(h) 
    set(h(j),'color',[1-(j-1)/(length(h)-1) 0 (j-1)/(length(h)-1)]); 
end