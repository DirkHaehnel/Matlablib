clear all
close all

if 0
    n = 1000;
    a = 1e-1;
    b = 2e-1;
    
    r = [0 0 0];
    v = [0 0 1];
    for j=2:n
        v = v + a*randn(1,3); 
        v = v/sqrt(v*v');
        r(end+1,:) = r(end,:) + b*v;
    end

    phi = 0:pi/10:2*pi;

    dv = diff(r);

    n1(:,1) = -dv(:,2);
    n1(:,2) = dv(:,1);
    n1(:,3) = zeros(size(dv,1),1);
    n1 = n1./(sqrt(sum(n1.^2,2))*[1 1 1]);

    n2(:,1) = dv(:,2).*n1(:,3)-dv(:,3).*n1(:,2);
    n2(:,2) = dv(:,3).*n1(:,1)-dv(:,1).*n1(:,3);
    n2(:,3) = dv(:,1).*n1(:,2)-dv(:,2).*n1(:,1);
    n2 = n2./(sqrt(sum(n2.^2,2))*[1 1 1]);

    r = (r(1:end-1,:)+r(2:end,:))/2;

    row = ones(size(phi));
    surfl(r(:,1)*row+n1(:,1)*cos(phi)+n2(:,1)*sin(phi),r(:,2)*row+n1(:,2)*cos(phi)+n2(:,2)*sin(phi),r(:,3)*row+n1(:,3)*cos(phi)+n2(:,3)*sin(phi));%ones(size(r,1)-1,size(phi)-1));
    shading flat
    axis image; axis off

end

if 0
    z = (0:0.01:0.8)';
    alpha = [0.3*cos(2*pi*z) 0.3*sin(2*pi*z) z];
    z = (0:0.01:0.6)'; phi = (0:pi/50:pi)';
    beta = [0*z 0*z z]; 
    beta  = [beta; [0.5*(1-cos(phi)) 0*phi z(end)+0.5*sin(phi)]];
    beta  = [beta; [ones(size(z)) 0*z z(end)-z]];
    
    n = 250;
	a = 1e-1;
    b = 1e-1;
    r = alpha;
    v = r(end,:)-r(end-1,:);
    v = v/sqrt(v*v');
    for j=2:50
        v = v + a*randn(1,3); 
        v = v/sqrt(v*v');
        r(end+1,:) = r(end,:) + b*v;
    end
	v = r(end,:)-r(end-1,:);
    v = v/sqrt(v*v');
    phi = angle(v(1)+i*v(2));
    theta = angle(v(3)+i*sqrt(sum(v(1:2).^2)));
    M = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)]*[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    r = [r; ones(size(beta,1),1)*r(end,:)+(M*beta')'];
    v = r(end,:)-r(end-1,:);
    v = v/sqrt(v*v');
    for j=2:50
        v = v + a*randn(1,3); 
        v = v/sqrt(v*v');
        r(end+1,:) = r(end,:) + b*v;
    end
	v = r(end,:)-r(end-1,:);
    v = v/sqrt(v*v');
    phi = angle(v(1)+i*v(2));
    theta = angle(v(3)+i*sqrt(sum(v(1:2).^2)));
    M = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)]*[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]*[1 0 0; 0 0 1; 0 -1 0];
    r = [r; ones(size(alpha,1),1)*(r(end,:)-(M*alpha(1,:)')')+(M*alpha')'];
	v = r(end,:)-r(end-1,:);
    v = v/sqrt(v*v');
    phi = angle(v(1)+i*v(2));
    theta = angle(v(3)+i*sqrt(sum(v(1:2).^2)));
    M = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)]*[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    r = [r; ones(size(beta,1),1)*r(end,:)+(M*beta')'];
    v = r(end,:)-r(end-1,:);
    v = v/sqrt(v*v');
    for j=2:50
        v = v + a*randn(1,3); 
        v = v/sqrt(v*v');
        r(end+1,:) = r(end,:) + b*v;
    end
	v = r(end,:)-r(end-1,:);
    v = v/sqrt(v*v');
    phi = angle(v(1)+i*v(2));
    theta = angle(v(3)+i*sqrt(sum(v(1:2).^2)));
    M = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)]*[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]*[1 0 0; 0 0 1; 0 -1 0];
    r = [r; ones(size(alpha,1),1)*(r(end,:)-(M*alpha(1,:)')')+(M*alpha')'];
    
    %plot3(r(:,1),r(:,2),r(:,3))
    
    phi = 0:pi/20:2*pi;

    dv = diff(r);

    n1(:,1) = -dv(:,2);
    n1(:,2) = dv(:,1);
    n1(:,3) = zeros(size(dv,1),1);
    n1 = 0.1*n1./(sqrt(sum(n1.^2,2))*[1 1 1]);

    n2(:,1) = dv(:,2).*n1(:,3)-dv(:,3).*n1(:,2);
    n2(:,2) = dv(:,3).*n1(:,1)-dv(:,1).*n1(:,3);
    n2(:,3) = dv(:,1).*n1(:,2)-dv(:,2).*n1(:,1);
    n2 = 0.1*n2./(sqrt(sum(n2.^2,2))*[1 1 1]);

    r = (r(1:end-1,:)+r(2:end,:))/2;

    row = ones(size(phi));
    surfl(r(:,1)*row+n1(:,1)*cos(phi)+n2(:,1)*sin(phi),r(:,2)*row+n1(:,2)*cos(phi)+n2(:,2)*sin(phi),r(:,3)*row+n1(:,3)*cos(phi)+n2(:,3)*sin(phi));%ones(size(r,1)-1,size(phi)-1));
    shading flat
    axis image; axis off

end

if 0
    z = (0:0.01:2.5)';
    alpha = [0.3*cos(2*pi*z) 0.3*sin(2*pi*z) z];
    z = (0:0.01:2)'; phi = (0:pi/20:pi)';
    beta = [0*z 0*z z]; 
    beta  = [beta; [0.5*(1-cos(phi)) 0*phi z(end)+0.5*sin(phi)]];
    beta  = [beta; [ones(size(z)) 0*z z(end)-z]];
    
    n = 50;
	a = 1e-1;
    b = 1e-2;
    r = alpha;
    v = r(end,:)-r(end-1,:);
    v = v/sqrt(v*v');
    for j=2:50
        v = v + a*randn(1,3); 
        v = v/sqrt(v*v');
        r(end+1,:) = r(end,:) + b*v;
    end
	v = r(end,:)-r(end-1,:);
    v = v/sqrt(v*v');
    phi = angle(v(1)+i*v(2));
    theta = angle(v(3)+i*sqrt(sum(v(1:2).^2)));
    M = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)]*[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    r = [r; ones(size(beta,1),1)*r(end,:)+(M*beta')'];
    v = r(end,:)-r(end-1,:);
    v = v/sqrt(v*v');
    for j=2:50
        v = v + a*randn(1,3); 
        v = v/sqrt(v*v');
        r(end+1,:) = r(end,:) + b*v;
    end
	v = r(end,:)-r(end-1,:);
    v = v/sqrt(v*v');
    phi = angle(v(1)+i*v(2));
    theta = angle(v(3)+i*sqrt(sum(v(1:2).^2)));
    M = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)]*[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]*[1 0 0; 0 0 1; 0 -1 0];
    r = [r; ones(size(alpha,1),1)*(r(end,:)-(M*alpha(1,:)')')+(M*alpha')'];
    
    phi = 0:pi/20:2*pi;

    dv = diff(r);

    n1(:,1) = -dv(:,2);
    n1(:,2) = dv(:,1);
    n1(:,3) = zeros(size(dv,1),1);
    n1 = 0.1*n1./(sqrt(sum(n1.^2,2))*[1 1 1]);

    n2(:,1) = dv(:,2).*n1(:,3)-dv(:,3).*n1(:,2);
    n2(:,2) = dv(:,3).*n1(:,1)-dv(:,1).*n1(:,3);
    n2(:,3) = dv(:,1).*n1(:,2)-dv(:,2).*n1(:,1);
    n2 = 0.1*n2./(sqrt(sum(n2.^2,2))*[1 1 1]);

    r = (r(1:end-1,:)+r(2:end,:))/2;

    row = ones(size(phi));
    surfl(r(:,1)*row+n1(:,1)*cos(phi)+n2(:,1)*sin(phi),r(:,2)*row+n1(:,2)*cos(phi)+n2(:,2)*sin(phi),r(:,3)*row+n1(:,3)*cos(phi)+n2(:,3)*sin(phi));%ones(size(r,1)-1,size(phi)-1));
    shading flat
    axis image; axis off

end

if 1
    dx = 0.001;
    x = 0:dx:1;
    z1 = 0.1*exp(-(x-0.4).^2/0.2^2);
    z2 = 0.9*exp(-(x-0.8).^2/0.1^2);
    tmp = sum(z1+z2)*dx;
    z1 = z1/tmp; z2 = z2/tmp;
    xm = sum(x.*z1+x.*z2)/sum(z1+z2);
    plot(x,z1+z2,'o',x,z1,x,z2,[xm xm],[0 5],'--');
    set(gca,'ytick',0:5,'fontsize',20)
    xlabel('FRET efficiency value','fontsize',24)
    ylabel('probability distribution','fontsize',24)
end

