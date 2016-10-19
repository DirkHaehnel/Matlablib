close;
n0 = 1.33;
n1 = 1.5;

if 0 % finding the refractive index
    theta=0:pi/5000:pi/2;
    xv = 0.63:0.001:0.68;
    yv = 2.6:0.01:3;
    ind = zeros(length(xv),length(yv));
    y = -0.5;
    for jx = 1:length(xv)
        x = xv(jx);
        for jy = 1:length(yv)
            y = yv(jy);
            [rp,rs] = DielMirror(theta,n0,n1,[n1 (n1+y)],[x*pi/n1 x*pi/(n1+y)],100);
            plot(theta,abs(rp),theta,abs(rs));
            axis([0 pi/2 0 1]);
            drawnow;
            if all(abs(rp)>1-1e-10 & abs(rs)>1-1e-10)
                ind(jx,jy)=1;
            end
        end
    end
    pcolor(n1+yv,xv,ind)
    xlabel('refractive index \itn\rm_2')
    ylabel('layer thickness (\lambda/\itn\rm)')
end

if 0 % calculating the bandgap
    theta=0:pi/5000:pi/2;
    lamem = 670;
    lambda = 550:0.2:750;
    x = 0.65;
    y = 3;
    for j=1:length(lambda)
        lam = lambda(j);
        [rp(j,:),rs(j,:)] = DielMirror(theta,n0,n1,[n1 (n1+y)],lam/lamem*[x*pi/n1 x*pi/(n1+y)],100);
    end
    ind = 1:length(lambda);
    ind = ind(all(abs(rp)>1-1e-4,2));
    pcolor(theta/pi*180,lambda,abs(rp));
    shading interp
	hold on
    line(theta/pi*180, lambda(ind(1))*ones(1,length(theta)), 'color', 'y', 'linewidth', 1)
    line(theta/pi*180, lambda(ind(end))*ones(1,length(theta)), 'color', 'y', 'linewidth', 1)
    hold off
    xlabel('incidence angle (°)');
    ylabel('wavelength (nm)');
    colorbar
    lammin = lambda(ind(1));
    lammax = lambda(ind(end));
    % print -dpng -r300 -zbuffer OmniRp
    % lammin = 651.4 nm
    % lammax = 687 nm
    if 0
        ind = 1:length(lambda);
        ind = ind(all(abs(rs)>1-1e-4,2));
        pcolor(theta/pi*180,lambda,abs(rs));
        shading interp
        hold on
        line(theta/pi*180, lambda(ind(1))*ones(1,length(theta)), 'color', 'y', 'linewidth', 1)
        line(theta/pi*180, lambda(ind(end))*ones(1,length(theta)), 'color', 'y', 'linewidth', 1)
        hold off
        xlabel('incidence angle (°)');
        ylabel('wavelength (nm)');
        colorbar
        % print -dpng -r300 -zbuffer OmniRs
        % lammax = 687 nm
    end
end

if 0 % calculating optimal excitation
    % lammin = 652 nm
    % lammax = 687 nm
    lamem = 670;
    lamex = 640;
    x = 0.65;
    y = 3;
    zv = 0:0.001:pi/5;
    clear rp
    for j=1:length(zv)
        rp(j) = Fresnel(n0,[n0 n1+y repmat([(n1+y) n1],1,50) n1],[zv(j) repmat(lamex/lamem*[x*pi/(n1+y) x*pi/n1],1,50)]);    
    end
    plot(zv/2/pi,abs(1-rp).^2)
    xlabel('high ref. index layer thickness (\lambda)')
    ylabel('rel. excitation intensity');
    zmax = zv(abs(1-rp).^2==max(abs(1-rp).^2));
    % zmax = 0.214;
end

if 0 % calcuating CEF (include distance dependence!) repeat with gold layer
    lammin = 652;
    lammax = 687;
    zmax = 0.214;
    dtheta = pi/1e5;
    theta = (asin(n0/n1)+dtheta/2):dtheta:pi/2;
    theta1 = dtheta/2:dtheta:pi/2;
    lamem = 670;
    lambda = lammin:0.02:lammax;
    x = 0.65;
    y = 3;
    zv = 0:pi/10:5*pi;
    ind1 = theta<asin(1.4/n1); % 1.4 NA objective
    ind2 = theta<asin(1.45/n1); % 1.45 NA objective    
    clear vres pres vres1 pres1 vres2 pres2
    for j=1:length(lambda)
        lam = lambda(j);
        for k=1:length(zv)
            z = zv(k);
            [v,pc,ps] = DipoleL(theta,z,[n1 repmat([(n1+y) n1],1,50) n1+y],n0,n0,[repmat(lam/lamem*[x*pi/(n1+y) x*pi/n1],1,50) zmax],max(z),[]);
            [v1,pc1,ps1] = DipoleL(theta1,z,n0,n0,[n1+y repmat([n1 (n1+y)],1,50) n1],[],max(z),[zmax repmat(lam/lamem*[x*pi/n1 x*pi/(n1+y)],1,50)]);
            [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu,qv,qp] = LifetimeL(z,[n1 repmat([(n1+y) n1],1,50) n1+y],n0,n0,[repmat(lam/lamem*[x*pi/(n1+y) x*pi/n1],1,50) zmax],max(z),[]);
            % fullangle:    
            vres(j,k) = dtheta*sin(theta)*abs(v).^2./(qvu+qvd+qv);
            pres(j,k) = 0.5*dtheta*(sin(theta)*(abs(pc).^2+abs(ps).^2))./(qpu+qpd+qp);
            
            % 1.4 NA:
            vres1(j,k) = dtheta*sin(theta(ind1))*abs(v(ind1)).^2./(qvu+qvd+qv);
            pres1(j,k) = 0.5*dtheta*(sin(theta(ind1))*(abs(pc(ind1)).^2+abs(ps(ind1)).^2))./(qpu+qpd+qp);
            
            % 1.45 NA:
            vres2(j,k) = dtheta*sin(theta(ind2))*abs(v(ind2)).^2./(qvu+qvd+qv);
            pres2(j,k) = 0.5*dtheta*(sin(theta(ind2))*(abs(pc(ind2)).^2+abs(ps(ind2)).^2))./(qpu+qpd+qp);
            j
        end
    end
end

if 0 % calculating angular distribution of radiation
    zmax = 0.214;
    theta = (asin(n0/n1)+pi/2e4):pi/1e5:pi/2;
    theta1 = pi/2e5:pi/1e5:pi/2;
    lamem = 670;
    x = 0.65;
    y = 3;
    z = 0;
    [v,pc,ps] = DipoleL(theta,z,[n1 repmat([(n1+y) n1],1,50) n1+y],n0,n0,[repmat(lam/lamem*[x*pi/(n1+y) x*pi/n1],1,50) zmax],max(z),[]);
    [v1,pc1,ps1] = DipoleL(theta1,z,n0,n0,[n1+y repmat([n1 (n1+y)],1,50) n1],[],max(z),[zmax repmat(lam/lamem*[x*pi/n1 x*pi/(n1+y)],1,50)]);
    plot(sin(theta).*abs(v)'.^2,-cos(theta).*abs(v)'.^2,-sin(theta).*abs(v)'.^2,-cos(theta).*abs(v)'.^2,'r',sin(theta1).*abs(v1)'.^2,cos(theta1).*abs(v1)'.^2,'b',-sin(theta1).*abs(v1)'.^2,cos(theta1).*abs(v1)'.^2,'b');
    plot(0.5*sin(theta).*(abs(pc).^2+abs(ps).^2)',-0.5*cos(theta).*(abs(pc).^2+abs(ps).^2)',-0.5*sin(theta).*(abs(pc).^2+abs(ps).^2)',-0.5*cos(theta).*(abs(pc).^2+abs(ps).^2)','r',0.5*sin(theta1).*(abs(pc1).^2+abs(ps1).^2)',0.5*cos(theta1).*(abs(pc1).^2+abs(ps1).^2)','b',-0.5*sin(theta1).*(abs(pc1).^2+abs(ps1).^2)',0.5*cos(theta1).*(abs(pc1).^2+abs(ps1).^2)','b')    
end

