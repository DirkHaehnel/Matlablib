close all

[x0,y0,z0] = sphere(30);
[x1,y1,z1] = sphere(50);

N = 200;
a = 5;
b = 25;
par = 3*[a a b];

for j=1:N
    pos(j,:) = 2*(rand(1,3)-0.5).*par;
end

res = [];
for k=1:200
    clf
    pos = pos + randn(size(pos));
    for j=1:N
        for jj = 1:3
            if pos(j,jj)>par(jj) pos(j,jj) = 2*par(jj) - pos(j,jj); end
            if pos(j,jj)<-par(jj) pos(j,jj) = -pos(j,jj) - 2*par(jj); end
        end
        col(j) = (sum(pos(j,1:2).^2)/a^2+pos(j,3)^2/b^2)<1;
    end
    subplot(121)
    hold on
    for j=1:N
        surf(x0 + pos(j,1),y0 + pos(j,2),z0 + pos(j,3), col(j)*ones(size(x0)))
    end
    surf(a*x1,a*y1,b*z1,0.7*ones(size(x1)));
    hold off
    axis image
    axis([-par(1) par(1) -par(2) par(2) -par(3) par(3)])
    box on
    set(gca,'xtick',[],'ytick',[],'ztick',[])
    caxis ([0 1])
    shading interp 
    camlight
    alpha(0.5)
    set(gcf,'color',[1 1 1])
    subplot(122)
    res = [res sum(col)];
    h = bar(0.5:k,res,1);
    set(h,'edgecolor','none');
    axis([0 200 0 10]);
    drawnow
    %eval(['print -dpng tmp' mint2str(k,3)]);
end
