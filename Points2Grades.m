function [g y n] = Points2Grades(x,x0,x1,y0,y1)

if nargin<4 || isempty(y0)
    y0 = 4;
end
if nargin<5 || isempty(y1)
    y1 = 1;
end

y = (y1-y0)/(x1-x0)*(x(:)-x0) + y0;
n = [1.0 1.3 1.7 2 2.3 2.7 3 3.3 3.7 4 4.3 4.7 5]';

ind = (n*ones(1,length(y)) - ones(length(n),1)*y').^2;
ind = ind==ones(size(ind,1),1)*min(ind);
g = 0*y;
for j=1:length(y)
    tmp = n(ind(:,j));
    g(j) = tmp(1);
end
g(g>4) = 5;
g(g<1) = 1;

z = mHist(g,n);
mHist(g,n)
hold on; 
for j=1:length(n) 
    text(n(j),z(j),int2str(z(j)),'horizontalalignment','center','verticalalignment','bottom'); 
end
hold off
set(gca,'XTick',n,'TickLength',[0 0])
xlabel('Note'); ylabel('Häufigkeit')
axis([0.7 5.3 0 max(z)+1.5])

