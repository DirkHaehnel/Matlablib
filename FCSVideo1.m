close all

[x0,y0,z0] = sphere(30);
[x1,y1,z1] = sphere(50);

N = 200;
a = 5;
b = 25;
par = 2*[b b b];
M = 200;

pos(1,:) = 2*(rand(1,3)-0.5).*par;
for j=2:N*M
    pos(j,:) = pos(j-1,:) + randn(1,3);
end
for j=1:N*M
    for jj = 1:3
        if pos(j,jj)>par(jj) pos(j,jj) = 2*par(jj) - pos(j,jj); end
        if pos(j,jj)<-par(jj) pos(j,jj) = -pos(j,jj) - 2*par(jj); end
    end
    col(j) = (sum(pos(j,1:2).^2)/a^2+pos(j,3)^2/b^2)<1;
end

for k=1:M
    clf
    hold on
    for j=1:N
        ind = (j-1)*M+k;        
        surf(x0 + pos(ind,1),y0 + pos(ind,2),z0 + pos(ind,3), col(ind)*ones(size(x0)))
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
    eval(['print -dpng tmp' mint2str(k,3)]);
end
