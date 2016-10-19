close all

[x0,y0,z0] = sphere(30);

N = 200;
b = 20;
par = 2*[b b b];
M = 200;

pos(1,:) = 2*(rand(1,3)-0.5).*par;
col(1) = 1/N;
for j=2:N*M
    pos(j,:) = pos(j-1,:) + randn(1,3);
    for jj = 1:3
        if pos(j,jj)>par(jj) pos(j,jj) = 2*par(jj) - pos(j,jj); end
        if pos(j,jj)<-par(jj) pos(j,jj) = -pos(j,jj) - 2*par(jj); end
    end
    col(j) = (mod(j-1,N)+1)/N;
end

for k=1:M
    clf
    hold on
    for j=1:N
        ind = (j-1)*M+k;        
        surf(x0 + pos(ind,1),y0 + pos(ind,2),z0 + pos(ind,3), col(ind)*ones(size(x0)))
    end
    hold off
    axis image
    axis([-par(1) par(1) -par(2) par(2) -par(3) par(3)])
    box on
    set(gca,'xtick',[],'ytick',[],'ztick',[])
    caxis ([0 5])
    shading interp 
    camlight
    alpha(0.5)
    set(gcf,'color',[1 1 1])
    drawnow
    eval(['print -dpng tmp' mint2str(k,3)]);
end
