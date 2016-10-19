% close
clear all

height = 10;
focus1 = 3;
lens1 = [150 height/sqrt(2)];
focus2 = 180;
lens2 = [height/sqrt(2) 500];
x0 = (-5:0.5:5)';
one = ones(size(x0));
p1 = height*one*[1,1]/sqrt(2) + [x0,one*10];
v1 = one*[0,-1];

for j=0:100
    theta = pi/4 + (j-50)/50*asin(5e-2/2/focus1);
    mirror = [cos(theta),sin(theta)];
    a(j+1) = theta-pi/4;

    phi = -pi/180*0;
    coverslide = [-cos(phi) sin(phi)];% -cos(phi)];

    lambda1 = (height-p1*mirror')./(v1*mirror');
    p2 = p1 + (lambda1*[1 1]).*v1;
    v2 = v1 - 2*(v1*mirror')*mirror;

    lambda2 = (one*lens1(1)-p2(:,1))./v2(:,1);
    p3 = p2 + (lambda2*[1 1]).*v2;
    v3(:,1) = one*focus1;
    v3(:,2) = v2(:,2)./v2(:,1)*focus1 - p3(:,2) + one*lens1(2);
    v3 = v3./(sqrt(sum(v3.^2,2))*[1 1]);
    
    lambda3 = (one*focus1)./v3(:,1);
    p4 = p3 + (lambda3*[1 1]).*v3;
    v4 = v3 - 2*(v3*coverslide')*coverslide;

    lambda4 = (one*focus1)./abs(v4(:,1));
    p5 = p4 + (lambda4*[1 1]).*v4;
    v5(:,1) = -one*focus1;
    v5(:,2) = -v4(:,2)./v4(:,1)*focus1 - p5(:,2) + one*lens1(2);
    v5 = v5./(sqrt(sum(v5.^2,2))*[1 1]);
    
    lambda5 = (height-p5*mirror')./(v5*mirror');
    p6 = p5 + (lambda5*[1 1]).*v5;
    v6 = v5 - 2*(v5*mirror')*mirror;

    lambda6 = (one*lens2(2)-p6(:,2))./v6(:,2);
    p7 = p6 + (lambda6*[1 1]).*v6;

    v7(:,2) = one*focus2;
    v7(:,1) = v6(:,1)./v6(:,2)*focus2 - p7(:,1) + one*lens2(1);
    v7 = v7./(sqrt(sum(v7.^2,2))*[1 1]);
    
    lambda7 = (one*focus2)./v7(:,2);
    p8 = p7 + (lambda7*[1 1]).*v7;

    % plot(p8(:,1)-one*height/sqrt(2), p8(:,2),'o');
    % plot(p8(:,1)-one*height/sqrt(2));
    % plot([p3(:,1) p4(:,1)]',[p3(:,2) p4(:,2)]');
    % plot([p7(:,1) p8(:,1)]',[p7(:,2) p8(:,2)]');
    % drawnow
    % figure(gcf)

    res(j+1,1) = mean(p4(:,2)-one*height/sqrt(2));
    res(j+1,2) = mean(p8(:,1)-one*height/sqrt(2));
end

plotyy(a,res(:,1),a,res(:,2))

return

plot([p1(:,1) p2(:,1)]',[p1(:,2) p2(:,2)]'); axis image; axis([0 153 -10 680]); drawnow
plot([p2(:,1) p3(:,1)]',[p2(:,2) p3(:,2)]'); axis image; axis([0 153 -10 680]); drawnow
plot([p3(:,1) p4(:,1)]',[p3(:,2) p4(:,2)]'); axis image; axis([0 153 -10 680]); drawnow
plot([p4(:,1) p5(:,1)]',[p4(:,2) p5(:,2)]'); axis image; axis([0 153 -10 680]); drawnow
plot([p5(:,1) p6(:,1)]',[p5(:,2) p6(:,2)]'); axis image; axis([0 153 -10 680]); drawnow
plot([p6(:,1) p7(:,1)]',[p6(:,2) p7(:,2)]'); axis image; axis([0 153 -10 680]); drawnow
plot([p7(:,1) p8(:,1)]',[p7(:,2) p8(:,2)]'); axis image; axis([0 153 -10 680]); drawnow
hold off
axis image
