clear all
close all

if 0 % Verteilung

    x=0:0.1:10;
    f=x.^2.*exp(-0.1*x.^2);f=f/sum(f)/mean(diff(x));

    plot(x,f)

    ind = x>4.8 & x<5.2; xx=x(ind); ff=f(ind);
    patch([xx(1:end-1); xx(2:end); xx(2:end); xx(1:end-1)],[0*ff(1:end-1); 0*ff(1:end-1); ff(2:end); ff(1:end-1)],'y','edgecolor','none')

    ind = x>1.7 & x<2.3; xx=x(ind); ff=f(ind);
    patch([xx(1:end-1); xx(2:end); xx(2:end); xx(1:end-1)],[0*ff(1:end-1); 0*ff(1:end-1); ff(2:end); ff(1:end-1)],'y','edgecolor','none')

    xlabel('\itx');
    ylabel('\itP\rm(\itx\rm)')
end

if 0 % Legendre Multipliers
    [x, y] = meshgrid(-1:0.01:1,-1:0.01:1);
    z0 = exp(-x.^2-y.^2);

    clear mov
    [xs,ys,zs] = sphere(20);
    cnt = 1;
    for lam = 0:0.01:0.45

        z = z0 - lam*(y+x+0.4);
        surf(x,y,z,'facealpha',0.5,'edgecolor','none');

        eps = 0.001;
        for j=-1:0.2:1
            line(x(abs(x-j)<eps),y(abs(x-j)<eps),z(abs(x-j)<eps),'color','k','linewidth',1);
            line(x(abs(y-j)<eps),y(abs(y-j)<eps),z(abs(y-j)<eps),'color','k','linewidth',1);
        end
        xx = -1:0.02:1;
        yy = -0.4-xx;
        xx = xx(abs(yy)<=1);
        yy = yy(abs(yy)<=1);

        line(xx,yy,exp(-xx.^2-yy.^2),'color','y');
        line(xx,yy,-1+0*xx,'color','y');        

        mx(cnt) = x(z==max(z(:))); my(cnt) = y(z==max(z(:))); mz(cnt) = z(z==max(z(:)));
        hold on
        surf(0.01*xs+mx(cnt),0.01*ys+my(cnt),0.01*zs+mz(cnt),0.6*ones(size(zs)),'edgecolor','none');

        line(mx,my,-1+0*mx,'color','r');                
        
        hold off

        axis image
        axis([-1 1 -1 1 -1 1]);
        set(gca,'linewidth',1,'xtick',[-1 0 1],'ytick',[-1 0 1],'ztick',[-1 0 1]);
        set(gcf,'color',[1 1 1])

        mov(cnt) = getframe; cnt = cnt+1;

    end

end

if 0 % Gibbs-Boltzmann
    clear all
    close all
    x= 1:100;
    n = 1000;
    beta = 0.05;
    n = n/sum(round(exp(-x*beta)));
    k = round(n*exp(-1*beta));
    plot([0 k+1],[1 1],1:k,ones(1,k),'o','linewidth',1); 
    hold on
    for j=2:length(x) 
        k = round(n*exp(-x(j)*beta));
        plot([0 k+1],[j j],1:k,j*ones(1,k),'o','linewidth',1); 
    end
    ylabel('\itE_j')
    axis([0 150 0 100])
end

if 0 % Maxwell
    clear all
    close all
    NatConst
    m = 4e-3/AvogadroConstant;
    
    x = 0:5e3;
    z1 = x.^2.*exp(-m/2*x.^2/BoltzmannConstant/(273.15+0)); z1 = z1/sum(z1)/mean(diff(x));
    z2 = x.^2.*exp(-m/2*x.^2/BoltzmannConstant/(273.15+100)); z2 = z2/sum(z2)/mean(diff(x));
    z3 = x.^2.*exp(-m/2*x.^2/BoltzmannConstant/(273.15+500)); z3 = z3/sum(z3)/mean(diff(x));
    plot(x,z1,x,z2,x,z3)
    
end

if 0 % Polymerkette End-2-end
    clear all
    close all
    NatConst
    a = 50;
    x = 1:1e3:3e5;
    plot(x,sqrt(x*a),'linewidth',3)
    set(gca,'fontsize',18)
%     xlabel('Konturlänge (nm)')
%     ylabel('\langle\itx\rm^2\rangle^{1/2} (nm)')
end

if 0 % Polymerkette Entropic Force
    clear all
    close all
    NatConst
    a = 50;
    L = 1e6*0.33;
    x = 1:1e3:L;
    f = 1:700;
    y = 1e-15*1e-9*f*a/BoltzmannConstant/(273.15+25);
    plot(x,3*1e15*1e9*BoltzmannConstant*(273.15+25)*x/L/a,L*sinh(y)./(2+cosh(y)),f,[L L],[0 max(f)],':','linewidth',3)
    set(gca,'fontsize',18)
end

if 0 % Kugel für Laplace
    [x,y,z]=sphere(50);
    surf(x,y,z,ones(size(x)))
    caxis([0 1.1])
    axis image
    axis off
    shading interp
    alpha(0.4)
    hold on
    pfeil([0.5 -0.5 0.5]/sqrt(3),[1 -1 1]/sqrt(3))
    pfeil([0.5 -0.5 0.5]/sqrt(3),[0 0 0])
    hold off
end

if 0 % Van-der-Waals-Gas
   a = 363.7*1e-3;
   b = 0.0427;
   V = (b+1e-5:0.001:0.5)';
   V = V/1e3; b= b/1e3;
   ax = [4000 9000]*1e3;
   T = 273.15 + [13 21 31 40 70];
   NatConst
   for j=1:length(T)
       p(:,j) = MolarGasConstant*T(j)./(V - b) - a./V.^2;
   end
   plot(V,p,'linewidth',3); set(gca,'fontsize',20); axis([min(V) max(V) ax]);
end

if 0 % chemical reaction 1
    NatConst
    x=0:0.001:1;
    temp = (273.15+25);
    plot(x,-5e3-MolarGasConstant*temp*log((2-2*x).^2.*(1-x)./x),x,0e3-MolarGasConstant*temp*log((2-2*x).^2.*(1-x)./x),...
        x,5e3-MolarGasConstant*temp*log((2-2*x).^2.*(1-x)./x),x,0*x,'--','linewidth',3); set(gca,'fontsize',20)
    grid
end

if 0 % chemical reaction 2
    NatConst
    x=0.0005:0.001:1;
    temp = (273.15+25);
    plot(x,(-5e3*x+MolarGasConstant*temp*(3*(1-x).*log(1-x)+x.*log(x) + 2*(1-log(2))*x))/1e3,...
        x,-5*x,'r:',...
        x,(0e3*x+MolarGasConstant*temp*(3*(1-x).*log(1-x)+x.*log(x) + 2*(1-log(2))*x))/1e3,...
        x,(5e3*x+MolarGasConstant*temp*(3*(1-x).*log(1-x)+x.*log(x) + 2*(1-log(2))*x))/1e3,...
        x,5*x,'g:',...
        x,0*x,'b:','linewidth',3); set(gca,'fontsize',20)
    plot(x,(-5e3*x+MolarGasConstant*temp*((2-2*x).*log(2-2*x)+(1-x).*log(1-x)+x.*log(x) + 2*(x-log(2))))/1e3,...
        x,-5*x,'r:',...
        x,(0e3*x+MolarGasConstant*temp*((2-2*x).*log(2-2*x)+(1-x).*log(1-x)+x.*log(x) + 2*(x-log(2))))/1e3,...
        x,(5e3*x+MolarGasConstant*temp*((2-2*x).*log(2-2*x)+(1-x).*log(1-x)+x.*log(x) + 2*(x-log(2))))/1e3,...
        x,5*x,'g:',...
        x,0*x,'b:','linewidth',3); set(gca,'fontsize',20)
    grid
end

if 0 % chemical reaction 1 gas
    NatConst
    x=0:0.001:1;
    temp = (273.15+25);
    plot(x,-5e3-MolarGasConstant*temp*log((2-2*x).^2.*(1-x)./x./(3-2*x).^2),...
        x,0e3-MolarGasConstant*temp*log((2-2*x).^2.*(1-x)./x./(3-2*x).^2),...
        x,5e3-MolarGasConstant*temp*log((2-2*x).^2.*(1-x)./x./(3-2*x).^2),...
        x,0*x,'--','linewidth',3); set(gca,'fontsize',20)
    grid
end

if 0 % chemical reaction 2 gas
    NatConst
    x=0.0005:0.001:1;
    temp = (273.15+25);
    plot(x,(-5e3*x+MolarGasConstant*temp*((2-2*x).*log(2-2*x)+(1-x).*log(1-x)+x.*log(x) - (3-2*x).*log(3-2*x) + (3*log(3)-2*log(2))))/1e3,...
        x,-5*x,'r:',...
        x,(0e3*x+MolarGasConstant*temp*((2-2*x).*log(2-2*x)+(1-x).*log(1-x)+x.*log(x) - (3-2*x).*log(3-2*x) + (3*log(3)-2*log(2))))/1e3,...
        x,(5e3*x+MolarGasConstant*temp*((2-2*x).*log(2-2*x)+(1-x).*log(1-x)+x.*log(x) - (3-2*x).*log(3-2*x) + (3*log(3)-2*log(2))))/1e3,...
        x,5*x,'g:',...
        x,0*x,'b:','linewidth',3); set(gca,'fontsize',20)
    grid
end

if 0 % HCN
    x = 10.^(-3:1e-5:0);
    Ks = 10^(-9.4);
    plot(x,(-Ks/2+sqrt((Ks/2)^2+Ks*x)),'linewidth',3); set(gca,'fontsize',20)
    figure;
    plot(x,(-Ks/2+sqrt((Ks/2)^2+Ks*x))./x,'b','linewidth',3); set(gca,'fontsize',20); set(gca,'YAxisLocation','right','XTick',[])
end

if 0 % Acetatpuffer
    x = -1:0.01:1;
    Ks = 10^(-4.75);
    cs = 1; cb = 1;
    z = -(cb+Ks+x)/2;
    z = z + sqrt(z.^2-cb*x+Ks*cs);
    plot(x,-log10(x+z),'linewidth',3); set(gca,'fontsize',20)
end

if 0 % electric field
   rho = 0.01*exp((0:0.01:1)*log(10));
   phi = 0:pi/300:2*pi;
   surf(rho'*cos(phi),rho'*sin(phi),1./rho'*ones(size(phi)));
   shading interp
   alpha(0.5)
   cameratoolbar
   hold on
   for j=1:12
       psi = 2*pi/12*(j-1);
       plot3(rho'*cos(psi),rho'*sin(psi),1./rho')
   end
   ax = axis;
   [x,y,z] = sphere(30);
   surf(0.01*x,0.01*y,0.01*z*diff(ax(5:6))/diff(ax(1:2))); shading interp
   hold off
end

if 0 % diffusion 1
   x = -10:0.01:10;
   t = 10.^(-2:1);
   j = 1;
   plot(x,1/sqrt(4*pi*t(j))*exp(-x.^2/4/t(j)));
   hold on
   for j = 2:length(t)
       plot(x,1/sqrt(4*pi*t(j))*exp(-x.^2/4/t(j)));
   end
   hold off
   colorize
end

if 0 % 2nd order reaction
    a0 = 1;
    b0 = 2;
    kp = 1;
    %km = 0.1;
    km = 0;
    k = kp/km;
    dc = b0 - a0;
    tmp = sqrt((dc+1/k)^2/4 + a0/k);
    lam1 = -(dc+1/k)/2 + tmp;
    lam2 = -(dc+1/k)/2 - tmp;
    t = 0:0.01:5;
    y = (a0-lam1)/(a0-lam2)*exp(-(lam1-lam2)*kp*t);
    y = (lam1-lam2*y)./(1-y);
    plot(t,y,t,b0-a0+y,t,a0-y,'linewidth',3); set(gca,'fontsize',20)
end

if 0 % 2nd order reaction special case
    a0 = 1;
    kp = 1;
    t = 0:0.01:5;
    y = a0./(1+a0*kp*t);
    plot(t,y,t,y,'r',t,a0-y,'linewidth',3); set(gca,'fontsize',20)
end

if 0 % 2nd order reaction comparison
    a0 = 1;
    b0 = 2;
    kp = 1;
    km = 0.1;
    k = kp/km;
    dc = b0 - a0;
    tmp = sqrt((dc+1/k)^2/4 + a0/k);
    lam1 = -(dc+1/k)/2 + tmp;
    lam2 = -(dc+1/k)/2 - tmp;
    t = 0:0.01:5;
    y1 = (a0-lam1)/(a0-lam2)*exp(-(lam1-lam2)*kp*t);
    y1 = (lam1-lam2*y1)./(1-y1);

    km = 0;
    k = kp/km;
    dc = b0 - a0;
    tmp = sqrt((dc+1/k)^2/4 + a0/k);
    lam1 = -(dc+1/k)/2 + tmp;
    lam2 = -(dc+1/k)/2 - tmp;
    y2 = (a0-lam1)/(a0-lam2)*exp(-(lam1-lam2)*kp*t);
    y2 = (lam1-lam2*y2)./(1-y2);

    y3 = a0./(1+a0*kp*t);
    
    plot(t,y1,t,y2,t,y3,'linewidth',3); set(gca,'fontsize',20)
end

if 0 % Michaelis menten
   k2 = 1;
   k = [0.5 1 2];
   cs = (0.01:0.01:10)';
   row = ones(size(k));
   col = ones(size(cs));
   plot(cs,(cs*row)./(cs*row+col*k),'linewidth',3); set(gca,'fontsize',20)
end

if 0 % Lineweaver-Burke 
   k2 = 0.2;
   k = [0.5 1 2];
   cs = (0.1:0.01:100)';
   row = ones(size(k));
   col = ones(size(cs));
   plot(1./cs,(cs*row+col*k)./(cs*row)/k2,'linewidth',3); set(gca,'fontsize',20)
end

if 0 % Eadie-Hofstee
   k2 = 0.2;
   k = [0.5 1 2];
   cs = (0.1:0.01:100)';
   row = ones(size(k));
   col = ones(size(cs));
   v = k2*(cs*row)./(cs*row+col*k);
   plot(v./(cs*row),v,'linewidth',3); set(gca,'fontsize',20)
end

if 0 % Hill coefficient
    k = 1;
    cs = (0.01:0.01:5)';
    plot(cs,cs./(cs+k),'linewidth',3); 
    hold on
    plot(cs,cs.^2./(cs.^2+k),'b','linewidth',3); 
    plot(cs,cs.^3./(cs.^3+k),'g','linewidth',3); 
    set(gca,'fontsize',20)
    hold off
end

if 0 % particle interference
    len = 200;
    [x,y] = meshgrid(0:0.2:len,-len:0.2:len);
    z = 20*pi;
    %     feld1 = exp(-(y+z/2).^2./50^2)./(1e2+x.^2+(y+z/2).^2);
    %     feld2 = exp(-(y-z/2).^2./50^2)./(1e2+x.^2+(y-z/2).^2);
    feld1 = 1./(1e2+x.^2+(y+z/2).^2);
    feld2 = 1./(1e2+x.^2+(y-z/2).^2);
    
    subplot(131); mim(log(feld1))
    subplot(132); mim(log(feld2))
    subplot(133); mim(log(feld1+feld2))
    figure
    plot(y(:,1), feld1(:,end)/max(feld1(:,end)), ...
        y(:,1), feld2(:,end)/max(feld1(:,end)), ...
        y(:,1), (feld1(:,end) + feld2(:,end))/max(feld1(:,end)),'linewidth',3); 
    set(gca,'fontsize',20)
end

if 0 % wave interference
    len = 200;
    [x,y] = meshgrid(0:0.2:len,-len:0.2:len);
    z = 20*pi;
    feld1 = exp(i*sqrt(x.^2+(y+z/2).^2))./(1e2+x.^2+(y+z/2).^2);
    feld2 = exp(i*sqrt(x.^2+(y-z/2).^2))./(1e2+x.^2+(y-z/2).^2);
    subplot(131); mim(log(abs(feld1)))
    subplot(132); mim(log(abs(feld2)))
    subplot(133); mim(log(abs(feld1+feld2)))

    figure
    plot(y(:,1), abs(feld1(:,end))/max(abs(feld1(:,end))), ...
        y(:,1), abs(feld2(:,end))/max(abs(feld1(:,end))), ...
        y(:,1), abs(feld1(:,end) + feld2(:,end))/max(abs(feld1(:,end))),'linewidth',3); 
    set(gca,'fontsize',20)
end

if 0 % water wave interference
    len = 200;
    [x,y] = meshgrid(0:0.2:len,-len:0.2:len);
    z = 20*pi;
    feld1 = exp(i*sqrt(x.^2+(y+z/2).^2));
    feld2 = exp(i*sqrt(x.^2+(y-z/2).^2));
    subplot(131); mim(real(feld1))
    subplot(132); mim(real(feld2))
    subplot(133); mim(real(feld1+feld2))
end

if 0 % Tunneleffekt
    NatConst
    v = 1e3;
    k = ElectronMass*v/PlanckConstant*2*pi;
    kappa = k;
    L = 200e-9;
    b = (k^2 + kappa^2)*sinh(kappa*L)/(2*i*k*kappa*cosh(kappa*L)+ (k^2-kappa^2)*sinh(kappa*L));
    c = k*(k+i*kappa)*exp(kappa*L)/(2*i*k*kappa*cosh(kappa*L)+ (k^2-kappa^2)*sinh(kappa*L));    
    d = -2*k*(k-i*kappa)/(-(k-i*kappa)^2+exp(2*L*kappa)*(k+i*kappa)^2);
    t = 4*i*k*kappa*exp(L*(kappa-i*k))/(-(k-i*kappa)^2+exp(2*L*kappa)*(k+i*kappa)^2);
    x1 = (-1e3:10:0)*1e-9;
    x2 = 0:1e-8:L;
    x3 = L:1e-8:1e-6;
    plot(1e6*x1,real(exp(i*k*x1)),'c:',1e6*x1,real(b*exp(-i*k*x1)),'c:',1e6*x1,real(exp(i*k*x1)+b*exp(-i*k*x1)),'r',...
        1e6*x2,real(c*exp(-kappa*x2)+d*exp(kappa*x2)),'b',1e6*x2,real(c*exp(-kappa*x2)),'c:',1e6*x2,real(d*exp(kappa*x2)),'c:',...
        1e6*x3,real(t*exp(i*k*x3)),'r','linewidth',3)    
    set(gca,'fontsize',20)
    axis tight
    figure
    plot(1e6*x1,abs(exp(i*k*x1)).^2,'c:',1e6*x1,abs(b*exp(-i*k*x1)).^2,'c:',1e6*x1,abs(exp(i*k*x1)+b*exp(-i*k*x1)).^2,'r',...
        1e6*x2,abs(c*exp(-kappa*x2)+d*exp(kappa*x2)).^2,'b',1e6*x2,abs(c*exp(-kappa*x2)).^2,'c:',1e6*x2,abs(d*exp(kappa*x2)).^2,'c:',...
        1e6*x3,abs(t*exp(i*k*x3)).^2,'r','linewidth',3)    
    set(gca,'fontsize',20)
    axis([-1 1 0 4])
end

if 0 % Einzelspalt
    NatConst
    v = 1e3;
    dv = 1e3; w = dv/2*ElectronMass*2*pi/PlanckConstant;
    y = (-5e3:10:5e3)*1e-9;
    x = (0:10:1e4)*1e-9;
    [x,y] = meshgrid(x,y);
    k = ElectronMass*v/PlanckConstant*2*pi;
    psi = exp(i*k*x-y.^2./(4/w^2+2*i*x/k));
    pcolor(1e6*x,1e6*y,real(psi)); shading interp; axis image; 
    figure
    pcolor(1e6*x,1e6*y,abs(psi).^2); shading interp; axis image
end

if 0 % Doppelspalt
    NatConst
    v = 1e3;
    dv = 1e3; w = dv/2*ElectronMass*2*pi/PlanckConstant;
    delta = 2e-6;
    y = (-5e3:10:5e3)*1e-9;
    x = (0:10:1e4)*1e-9;
    [x,y] = meshgrid(x,y);
    k = ElectronMass*v/PlanckConstant*2*pi;
    psi = exp(i*k*x-(y-delta).^2./(4/w^2+2*i*x/k)) + exp(i*k*x-(y+delta).^2./(4/w^2+2*i*x/k));
    pcolor(1e6*x,1e6*y,real(psi)); shading interp; axis image; 
    figure
    pcolor(1e6*x,1e6*y,abs(psi).^2); shading interp; axis image
end

if 0 % Fourier
    NatConst
    sig = 1e-7;
    %k = 10.^(8:-0.025:5);    
    k = 10.^(5:0.025:8);    
    dk = gradient(k);
    x = (-1e3:1:1e3)*1e-9;
    tst = 0;
    for j = 1:length(k)
        tst = tst + abs(dk(j))*exp(i*x*k(j) - k(j)^2*sig^2);
        plot(1e6*x,real(exp(i*x*k(j))),1e6*x,abs(tst).^2/max(abs(tst).^2))
        axis([-1 1 -1 1])
        set(gca,'xtick',[],'ytick',[])
        mov(j) = getframe;
    end
end

if 0 % Wellenpaket
    NatConst
    v = 0e4;
	sig = 1e-7;
    x = (-1e4:10:1e4)*1e-9;
    t = 0:1e-11:5e-10;
    clear mov
    for j=1:length(t)
        %psi = 1/(2*pi*sig^2*(1+PlanckConstant*t(j)/2/ElectronMass/sig^2)^2)^(1/4)*exp(-(x-v*t(j)).^2/4/(sig^2+i*PlanckConstant*t(j)/2/ElectronMass));
        psi = exp(-(x-v*t(j)+v*max(t)/2).^2/4/(sig^2+i*PlanckConstant*t(j)/2/ElectronMass));
        plot(1e6*x,real(psi),1e6*x,abs(psi).^2,'linewidth',3)    
        set(gca,'fontsize',20)
        axis([-10 10 -1 1])
        set(gca,'xtick',[],'ytick',[])
        mov(j) = getframe;
    end
end

if 0 % Potentialwall
    NatConst
    L = 200e-9;
    ne = 0:2;
    nu = 1:3;
    ke = (ne+0.5)*2*pi/L;
    ku = nu*2*pi/L;
    x = (-L/2:1e-9:L/2);
    ee = (PlanckConstant*ke).^2/2/ElectronMass;
    eu = (PlanckConstant*ku).^2/2/ElectronMass;    
    
    subplot(121)
    plot(1e9*x,2e22*(ee'*ones(size(x))) + cos(ke'*x),'b',1e9*x,2e22*(eu'*ones(size(x))) + sin(ku'*x),'r','linewidth',3)
    hold on
    plot(1e9*x,2e22*(ee'*ones(size(x))),'k',1e9*x,2e22*(eu'*ones(size(x))),'k','linewidth',1)
    hold off
    set(gca,'fontsize',20)
    subplot(122)
    plot(1e9*x,2e22*(ee'*ones(size(x))) + abs(cos(ke'*x)).^2,'b',1e9*x,2e22*(eu'*ones(size(x))) + abs(sin(ku'*x)).^2,'r','linewidth',3)
    hold on
    plot(1e9*x,2e22*(ee'*ones(size(x))),'k',1e9*x,2e22*(eu'*ones(size(x))),'k','linewidth',1)
    hold off
    set(gca,'fontsize',20)
end

if 0 % Harmonischer Oscillator Grundzustand
    NatConst
    L = 0.7e-11;
    ne = 0:5;
    m = 12*ProtonMass;
    om = 2*pi*SpeedOfLight/1e-5;
    x = 5*(-L/2:1e-13:L/2);
    eps = m*om/PlanckConstant*2*pi;
    
    plot(1e9*x,(m*om.^2)/2*x.^2,'b',1e9*x,(PlanckConstant/2/pi*om/2)*(1+0.5*exp(-eps*x.^2/2)),'r','linewidth',3)
    hold on
    plot(1e9*x,(PlanckConstant/2/pi*om/2)*ones(size(x)),'k','linewidth',1)
    hold off
    axis([-5/2*L*1e9 5/2*L*1e9 0 2e-20])
    set(gca,'fontsize',20)
end

if 0 % Harmonischer Oscillator 
    NatConst
    L = 2e-11;
    n = 5;
    m = 12*ProtonMass;
    om = 2*pi*SpeedOfLight/1e-5;
    x = 5*(-L/2:1e-13:L/2);
    eps = m*om/PlanckConstant*2*pi;

    subplot(121)    
    plot(1e9*x,(m*om.^2)/2*x.^2,'b',1e9*x,(PlanckConstant/2/pi*om/2)*(1+0.5*exp(-eps*x.^2/2)),'r','linewidth',3)
    hold on
    plot(1e9*x,(PlanckConstant/2/pi*om/2)*ones(size(x)),'k','linewidth',1)
    for j=1:n
        tst = exp(-eps*x.^2/2).*polyval(HermitePoly(j),sqrt(eps)*x);
        tst = tst/max(tst)*(PlanckConstant/2/pi*om/2)*0.5;
        plot(1e9*x,(PlanckConstant/2/pi*om*(j+0.5))+tst,'r','linewidth',3)
        hold on
        plot(1e9*x,(PlanckConstant/2/pi*om*(j+0.5))*ones(size(x)),'k','linewidth',1)
    end
    hold off
    axis([-0.025 0.025 0 1.5e-19])
    set(gca,'fontsize',20)
    
    subplot(122)    
    plot(1e9*x,(m*om.^2)/2*x.^2,'b',1e9*x,(PlanckConstant/2/pi*om/2)*(1+0.5*abs(exp(-eps*x.^2/2)).^2),'r','linewidth',3)
    hold on
    plot(1e9*x,(PlanckConstant/2/pi*om/2)*ones(size(x)),'k','linewidth',1)
    for j=1:n
        tst = abs(exp(-eps*x.^2/2).*polyval(HermitePoly(j),sqrt(eps)*x)).^2;
        tst = tst/max(tst)*(PlanckConstant/2/pi*om/2)*0.5;
        plot(1e9*x,(PlanckConstant/2/pi*om*(j+0.5))+tst,'r','linewidth',3)
        hold on
        plot(1e9*x,(PlanckConstant/2/pi*om*(j+0.5))*ones(size(x)),'k','linewidth',1)
    end
    hold off
    axis([-0.025 0.025 0 1.5e-19])
    set(gca,'fontsize',20)
end

if 0 % Drehimpuls
    phi=0:pi/200:2*pi;
    surf(ones(2,1)*cos(phi),ones(2,1)*sin(phi),[1 -1]'*ones(size(phi))); 
    shading interp; 
    alpha(0.5); 
    line(cos(phi),sin(phi),0.2*cos(10*phi),'linewidth',3)
    line([0 0],[0 0],[-3 3],'color','b','linewidth',2)
    axis off
end

if 0 % Drehimpulsdrehung
    phi=0:pi/200:2*pi;
    subplot(131)
    surf(ones(2,1)*cos(phi),ones(2,1)*sin(phi),0.5*[1 -1]'*ones(size(phi))); 
    shading interp; 
    alpha(0.5); 
    line(cos(phi),sin(phi),0.2*cos(10*phi),'linewidth',3)
    line([0 0],[0 0],[-1 1],'color','b','linewidth',2)
    axis off
    axis image

    subplot(132)
    surf(ones(2,1)*cos(phi),0.5*[1 -1]'*ones(size(phi)),ones(2,1)*sin(phi)); 
    shading interp; 
    alpha(0.5); 
    line(cos(phi),0.2*cos(10*phi),sin(phi),'linewidth',3)
    line([0 0],[-1 1],[0 0],'color','b','linewidth',2)
    axis off
    axis image
    
    subplot(133)
    surf(ones(2,1)*cos(phi),ones(2,1)*sin(phi),0.5*[1 -1]'*ones(size(phi))); 
    shading interp; 
    alpha(0.5); 
    line(cos(phi),sin(phi),-0.2*cos(10*phi),'linewidth',3)
    line([0 0],[0 0],[-1 1],'color','b','linewidth',2)
    axis off
    axis image    
end

if 0 % Hydrogen theta-phi
    theta = (0:pi/100:pi)';
    phi = 0:pi/100:2*pi;
    cnt = 1;
    
    [f theta phi] = ylm(0,0); theta = theta'; f = real(f);
    surf((sin(theta)*cos(phi)).*abs(f).^2,(sin(theta)*sin(phi)).*abs(f).^2,(cos(theta)*ones(size(phi))).*abs(f).^2,...
        permute(repmat([1 0 0]',[1 length(theta) length(phi)]),[2 3 1]));
    shading interp; alpha(0.5);
    a = permute(cat(3,zeros(3),0.15*eye(3)),[3,1,2]);
    line(a(:,:,1),a(:,:,2),a(:,:,3),'color',[0.5 0.3 1])
    view([137.5 30])
    axis off
    axis image
    eval(['print -dpng -r300 tmp' mint2str(cnt,2)]); cnt = cnt+1;

    for j=-1:1
        figure
        [f theta phi] = ylm(1,j); theta = theta'; f = real(f);
        surf((sin(theta)*cos(phi)).*abs(f).^2,(sin(theta)*sin(phi)).*abs(f).^2,(cos(theta)*ones(size(phi))).*abs(f).^2,f);
        shading interp; alpha(0.5);
        a = permute(cat(3,zeros(3),0.15*eye(3)),[3,1,2]);
        line(a(:,:,1),a(:,:,2),a(:,:,3),'color',[0.5 0.3 1])
        view([137.5 30])
        axis off
        axis image
        eval(['print -dpng -r300 tmp' mint2str(cnt,2)]); cnt = cnt+1;
    end

    for j=-2:2
        figure
        [f theta phi] = ylm(2,j); theta = theta'; f = real(f);
        surf((sin(theta)*cos(phi)).*abs(f).^2,(sin(theta)*sin(phi)).*abs(f).^2,(cos(theta)*ones(size(phi))).*abs(f).^2,f);
        shading interp; alpha(0.5);
        a = permute(cat(3,zeros(3),0.15*eye(3)),[3,1,2]);
        line(a(:,:,1),a(:,:,2),a(:,:,3),'color',[0.5 0.3 1])
        view([137.5 30])
        axis off
        axis image
        eval(['print -dpng -r300 tmp' mint2str(cnt,2)]); cnt = cnt+1;
    end

end

if 0 % Kugelkoordinaten

    theta = pi/4; phi = pi/3;
    
    [x,y,z] = sphere(100);
	surf(x,y,z,permute(repmat([1 0.8 0.8]',[1 size(x,1) size(x,2)]),[2 3 1]));
    shading interp
    axis image
    alpha(0.9);
    a = permute(cat(3,zeros(3),1.2*eye(3)),[3,1,2]);
    line(a(:,:,1),a(:,:,2),a(:,:,3),'color',[0.5 0.3 1])
    line(1.3*[0 sin(theta)*cos(phi)],1.3*[0 sin(theta)*sin(phi)],1.3*[0 cos(theta)],'color',[0.5 0.3 1]);
    line(1.3*[0 sin(theta)*cos(0)],1.3*[0 sin(theta)*sin(0)],1.3*[0 cos(theta)],'color',[0.5 0.3 1]);
    line(sin(0:pi/100:theta)*cos(phi),sin(0:pi/100:theta)*sin(phi),cos(0:pi/100:theta),'color',[0.5 0.3 1]);
    line(cos(0:pi/100:phi)*sin(theta),sin(0:pi/100:phi)*sin(theta),ones(size(0:pi/100:phi))*cos(theta),'color',[0.5 0.3 1]);
    line(sin(0:pi/100:pi),sin(0:pi/100:pi)*0,cos(0:pi/100:pi),'color',[1 0.5 0.5]);
    line(cos(0:pi/100:2*pi),sin(0:pi/100:2*pi),sin(0:pi/100:2*pi)*0,'color',[1 0.5 0.5]);

    view([137.5 30])
    axis off
    axis image
    camlight left

end

if 0 % Hydrogen radial
    N = 1;
    
    r = (0.005:0.01:12);
    
    subplot(323)
    L = 0;    
    tmp = exp(-2*r/N).*(2*r/N).^L.*polyval(AssociatedLaguerrePoly(N-L-1,2*L+1),2*r/N); tmp = tmp/max(tmp);
    plot(r,tmp)
    hold on
    for L=1:N-1
        tmp = exp(-2*r/N).*(2*r/N).^L.*polyval(AssociatedLaguerrePoly(N-L-1,2*L+1),2*r/N); tmp = tmp/max(tmp);
        plot(r,tmp)
    end
    hold off
    if N>1 colorize; end
    axis tight
    subplot(324)
    L = 0;    
    tmp = abs(exp(-2*r/N).*(2*r/N).^L.*polyval(AssociatedLaguerrePoly(N-L-1,2*L+1),2*r/N)).^2; tmp = tmp/max(tmp);
    plot(r,tmp)
    hold on
    for L=1:N-1
        tmp = abs(exp(-2*r/N).*(2*r/N).^L.*polyval(AssociatedLaguerrePoly(N-L-1,2*L+1),2*r/N)).^2; tmp = tmp/max(tmp);
        plot(r,tmp)
    end
    hold off
    if N>1 colorize; end
    axis tight
end

if 0 % Hydrogen energies
    NatConst
    m = 1/(1/ElectronMass + 1/ProtonMass);
    eps = -m*ElectronCharge^3/8/VacuumPermittivity^2/PlanckConstant^2;
    N = 10;
    
    for j=1:N
        for k=1:j
            line((k-1)*(1+0.5) + (0:1), eps/j^2*[1 1])
        end
    end
    set(gca,'xtick',[],'ytick',[],'yscale','log')
    axis([0 15 -20 -0.1])
end

if 0 % periodic table
    NatConst
    m = 1/(1/ElectronMass + 1/ProtonMass);
    eps = -m*ElectronCharge^3/8/VacuumPermittivity^2/PlanckConstant^2;
    N = 7;

    for j=1:N
        line((0:1), eps/j^2*[1 1])
    end
    set(gca,'xtick',[],'ytick',[],'yscale','log')
    axis([0 15 -20 -0.1])
end

if 0 % Benzen Energie
    phi = 0:pi/100:2*pi; phase = 0;
    plot(cos(phi),sin(phi),'-',cos(2*pi/6*(1:6)-phase),sin(2*pi/6*(1:6)-phase),'o',0,0,'.','linewidth',3)
    axis image
    axis off
end

if 0 % Benzen Wellenfunktion
    phi = 0:pi/100:2*pi; 
    delta = pi/30;
    for k = 1:6
        subplot(2,3,k)
        plot3(cos(phi),sin(phi),0*cos(phi));
        hold on
        plot3(cos(2*pi/6*(0:5)),sin(2*pi/6*(0:5)),zeros(1,6),'or','markersize',5);
        for j=1:6
            patch(cos(2*pi*j/6+[-delta -delta delta delta])',sin(2*pi*j/6+[-delta -delta delta delta])',...
                cos(k*2*pi*j/6)*[0 1 1 0]','b')
            alpha(0.2)
        end
        plot3(cos(phi),sin(phi),cos(k*phi),'b','linewidth',1);
        axis image
        axis off
    end
end

if 1 % Photon
    x = 0:pi/25:10*pi;
    plot3(x,0*x,0*x,'b');
    line(x,cos(x),-sin(x),'linewidth',3);
    axis image
    axis off
end