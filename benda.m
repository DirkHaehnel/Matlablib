clear all
if 1
    lam = 520;
    nv = [1.457 1.52 1.615];
    nspb = 1.45;
    nw = 1.333;
    d1 = 3/lam*2*pi;
    d = 3/lam*2*pi;
    d2 = [];
    z = 0;
    for k=1:3
        [lvd, lvu, lpd, lpu, qvd(k), qvu(k), qpd(k), qpu(k)] = LifetimeL(z,[nv(k) nspb],nspb,nw,d1,d,d2);
    end
    meas = [0.87 0.835 0.77];
    [lvd, lvu, lpd, lpu] = LifetimeL(d,nw,nspb,nw,[],2*d,[]);
    theta = (0:pi/200:pi/2)';
    tst = sum((((cos(theta).^2*(lvu+lvd)+sin(theta).^2*(lpu+lpd))*[1 1 1])./(cos(theta).^2*(qvu+qvd)+sin(theta).^2*(qpu+qpd))-ones(size(theta))*meas).^2,2);
    plot(theta/pi*180,((cos(theta).^2*(lvu+lvd)+sin(theta).^2*(lpu+lpd))*[1 1 1])./(cos(theta).^2*(qvu+qvd)+sin(theta).^2*(qpu+qpd)),...
        theta/pi*180,ones(size(theta))*meas,':');
    ax = axis;
    plot(theta/pi*180,((cos(theta).^2*(lvu+lvd)+sin(theta).^2*(lpu+lpd))*[1 1 1])./(cos(theta).^2*(qvu+qvd)+sin(theta).^2*(qpu+qpd)),...
        theta/pi*180,ones(size(theta))*meas,':',theta(tst==min(tst))/pi*180*[1 1],ax(3:4),':');
    xlabel('inclination angle [°]'); ylabel('rel. lifetime')    
    
end

if 0
    lam = 0.67;
    % n1 = [1.51, 0.3257 + 2.5792i];
    n1 = [1.52 2+0.03i]; n = 1.33; n2 = 1.33;
    d1 = 0.2/lam*2*pi; d = 0.2/lam*2*pi; d2 = [];
    z = (1:200)/1000/lam*2*pi;

    for j=1:length(z)
        [lvd(j), lvu(j), lpd(j), lpu(j), qvd(j), qvu(j), qpd(j), qpu(j)] = LifetimeL(z(j),n1,n,n2,d1,d,d2);
    end
    plot(1e3*z/2/pi*lam,4/3*n./(qvu+qvd),1e3*z/2/pi*lam,4/3*n./(qpu+qpd));
end

if 0
    lam = 0.67;
    n1 = [1.52 2+0.03i 1.33]; n = 1.4; n2 = 1.33;
    d1 = 0.2/lam*2*pi; d = 0.002/lam*2*pi; d2 = [];
    z = (1:200)/1000/lam*2*pi;

    for j=1:length(z)
        [lvd1(j), lvu1(j), lpd1(j), lpu1(j), qvd1(j), qvu1(j), qpd1(j), qpu1(j)] = LifetimeL(d/2,n1,n,n2,[d1 z(j)],d,d2);
    end
    plot(1e3*z/2/pi*lam,4/3*n./(qvu1+qvd1),1e3*z/2/pi*lam,4/3*n./(qpu1+qpd1));
end
