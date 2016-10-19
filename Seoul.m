% Seoul

phi=0:pi/50:3/2*pi;
theta=(0:pi/50:pi)';
rv = [1 1.2];

for k=1:10
    alpha = rand*2*pi;
    beta = rand*pi;
    for j=1:20
        r = rv(1);
        surf(r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta)*ones(size(phi)),0.2*ones(size(theta))*ones(size(phi)))
        hold on
        r = rv(2);
        surf(r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta)*ones(size(phi)),0.2*ones(size(theta))*ones(size(phi)))
        surf(sin(theta)*rv,zeros(length(theta),2),cos(theta)*rv,ones(length(theta),2)*0.2)
        surf(zeros(length(theta),2),-sin(theta)*rv,cos(theta)*rv,ones(length(theta),2)*0.2)
        [x,y] = meshgrid(-1.5:0.1:1.5,-1.5:0.1:1.5);
        delta = -2+j*0.2;
        z = -sin(beta)*x+cos(beta)*delta;
        x = cos(beta)*x+sin(beta)*delta;
        tmp = cos(alpha)*x+sin(alpha)*y;
        y = -sin(alpha)*x+cos(alpha)*y;
        x = tmp;
        surf(x,y,z,ones(size(x)))
        hold off
        axis image
        shading interp
        camlight
        axis off
        view([20 10])
        axis([-3 3 -3 3 -3 3])
        caxis([-0.2 1.2])
        mov(j+(k-1)*20)=getframe;
    end
end