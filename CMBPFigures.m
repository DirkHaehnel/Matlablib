close all

if 0
    
    x = [0 exp((0:1e2)*log(1e3)/100)];
    phi=0:pi/100:3/2*pi;
    psi=0:pi/100:2*pi;
    
    sig0=200;
    sig1=30;
    surf(x'*cos(phi),x'*sin(phi),exp(-x.^2/2/sig0^2)'*ones(size(phi)))
    hold on
    surf(x'*cos(phi),x'*sin(phi),exp(-x.^2/2/sig1^2)'*ones(size(phi)))
    plot(sig0*cos(psi),sig0*sin(psi),'b')
    plot(sig1*cos(psi),sig1*sin(psi),'r')
    hold off
    shading interp; alpha(0.5)
    axis off
    
end