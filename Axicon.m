w = 1000; 
z = 0; 
lambda = 1.064; 
Nmax = 600;
alpha = 1e-2;
dx = 1/7000;

x = (0.5:(1/dx))'*dx*60;
dx = diff(x(1:2));
rho = sqrt(x/2)*w*sqrt(1+(lambda*z/pi/w^2)^2);

pl = morthpl(3,Nmax,x);

tmp = exp(-i/2*alpha*rho).*pl(:,1).*exp(-x);
tmp = tmp.';

coef = dx*tmp*pl;

close 
for z=0:1e4:2e6
    rho = sqrt(x/2)*w*sqrt(1+(lambda*z/pi/w^2)^2);
    plot(rho,abs((pl*(coef.'.*exp(-i*(2*(0:Nmax)'+1)*atan(lambda*z/pi/w^2)))).*exp(-rho.^2/w^2/(1+(lambda*z/pi/w^2)^2))*2/w/sqrt(1+(lambda*z/pi/w^2)^2)).^2,rho,abs(exp(-x/2)*2/w).^2)
    %plot(sqrt(x/2),abs((pl*(coef.'.*exp(-i*(2*(0:Nmax)'+1)*atan(lambda*z/pi/w^2)))).*exp(-rho.^2/w^2/(1+(lambda*z/pi/w^2)^2))*2/w/sqrt(1+(lambda*z/pi/w^2)^2)).^2,sqrt(x/2),abs(exp(-x/2)*2/w).^2)
    title([num2str(z/1e4,3) ' cm'])
    %frame = getframe(gca);
    ax = axis;
    axis([ax(1) ax(2) 0 3e-5])
    drawnow
end



