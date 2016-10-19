function z = RotoDiffTop(theta,t,drot)

nn = 10;

dm = (drot(1)-drot(2))/2;
dp = (drot(1)+drot(2))/2;
dc = drot(3);

A = zeros(nn); 
for j=1:nn
    A(j,j) = 6*dp + dc*(2*j-2)^2;
    if j<nn
        A(j,j+1) = -3*dm;
    end
    if j>1 
        A(j,j-1) = -3*dm;
    else
        A(j,j+1) = -6*dm;
    end
end
[v, d] = eig(A);
vi = inv(v);

z = 1+2*sum(vi(:,1).*v(1,:)'.*exp(-diag(d)*t))*polyval(LegendrePoly(2),cos(theta));

return

close all; 
t=0:0.01:1; 
for j=1:length(t) 
    z0 = RotoDiffTop(theta,t(j),[1 1 1]);     
    z1 = RotoDiffTop(theta,t(j),[0.1 1 1]); 
    plot(z0.*cos(theta),z0.*sin(theta),'r',z0.*cos(theta),-z0.*sin(theta),'r',...
        z1.*cos(theta),z1.*sin(theta),'b',z1.*cos(theta),-z1.*sin(theta),'b'); 
    axis image; 
    axis([-3 3 -1.5 1.5]); 
    drawnow; 
end