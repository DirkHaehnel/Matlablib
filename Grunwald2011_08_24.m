clear all

nm = 1.33:0.0001:1.515; % sample ref. index
ng = [1.33 1.405 1.515]; % ref. index immersion medium of objective
NA = [1.2 1.3 1.49];%ng*sin(61.09/180*pi); % num. ap.
lamem = 0.58;
zp = 3/lamem*2*pi; % vector of considered distances between molecules and surface

if 1
    close all
    for j=1:numel(NA)
        % calculate max. collection angle
        theta_max = asin(NA(j)/ng(j));
        
        dtheta = pi/1e4;
        theta = dtheta/2:dtheta:pi;
        ind = theta<=theta_max;
        
        for k=1:numel(nm)
            % calculate angular dist. of emission into glass
            [v, pc, ps] = DipoleL(theta,zp,ng(j),nm(k),nm(k),[],zp,[]);
            
            % calculate angular dist. of emission into sample
            [v1,pc1,ps1] = DipoleL(theta,0,nm(k),nm(k),ng(j),[],zp,[]);
            
            % calculate collection efficiency
            cefv(k,j) = (sin(theta(ind))*abs(v(ind)).^2)/(sin(theta)*(abs(v).^2+abs(v1).^2));
            cefp(k,j) = (sin(theta(ind))*(abs(pc(ind)).^2+abs(ps(ind)).^2))/...
                (sin(theta)*(abs(pc).^2+abs(ps).^2+abs(pc1).^2+abs(ps1).^2));
            totalv(k,j) = (sin(theta)*abs(v).^2)/(sin(theta)*(abs(v).^2+abs(v1).^2));
            totalp(k,j) = (sin(theta)*(abs(pc).^2+abs(ps).^2))/...
                (sin(theta)*(abs(pc).^2+abs(ps).^2+abs(pc1).^2+abs(ps1).^2));
        end
    end
    
    plot(nm,(cefv+2*cefp)/3)
    xlabel('R.I.');
    ylabel('collection efficiency (a.u.)');
    legend({'NA = 1.2','NA = 1.3','NA = 1.49'},3)
    figure
    plot(nm,(totalv+2*totalp)/3)
    xlabel('R.I.');
    ylabel('total emission into glass (a.u.)');
    legend({'NA = 1.2','NA = 1.3','NA = 1.49'},3)
end

if 0
    close all
    zp = 3/lamem*2*pi;
    nmv = 1.33:0.005:1.45;
    for j=1:numel(nmv)
        nm = nmv(j);
        NA = 1.3;
        ng = 1.405;
        theta_max = asin(NA/ng);
        dtheta = pi/1e4;
        theta = (dtheta/2:dtheta:pi)';
        ind = theta<=theta_max;
        [v, pc, ps] = DipoleL(theta,zp,ng,nm,nm,[],zp,[]);
        [v1,pc1,ps1] = DipoleL(theta,0,nm,nm,ng,[],zp,[]);
        tst(j,:) = [sin(theta')*(abs(pc).^2+abs(ps).^2+abs(v).^2) sin(theta')*(abs(pc1).^2+abs(ps1).^2+abs(v1).^2)]*dtheta/2*3/2;
        nrmp = sqrt(sin(theta')*(abs(pc).^2+abs(ps).^2+abs(pc1).^2+abs(ps1).^2)*dtheta/2*3/2);
        nrmv = sqrt(sin(theta')*(abs(v).^2+abs(v1).^2)*dtheta*3/2);
        v = 0*v/nrmv; pc = pc/nrmp; ps = ps/nrmp;
        v1 = 0*v1/nrmv; pc1 = pc1/nrmp; ps1 = ps1/nrmp;
        plot(-sin(theta(ind)).*(abs(pc(ind)).^2+abs(ps(ind)).^2+abs(v(ind)).^2), -cos(theta(ind)).*(abs(pc(ind)).^2+abs(ps(ind)).^2+abs(v(ind)).^2),'r',...
            sin(theta(ind)).*(abs(pc(ind)).^2+abs(ps(ind)).^2+abs(v(ind)).^2), -cos(theta(ind)).*(abs(pc(ind)).^2+abs(ps(ind)).^2+abs(v(ind)).^2),'r',...
            -sin(theta(~ind)).*(abs(pc(~ind)).^2+abs(ps(~ind)).^2+abs(v(~ind)).^2), -cos(theta(~ind)).*(abs(pc(~ind)).^2+abs(ps(~ind)).^2+abs(v(~ind)).^2),'b',...
            sin(theta(~ind)).*(abs(pc(~ind)).^2+abs(ps(~ind)).^2+abs(v(~ind)).^2), -cos(theta(~ind)).*(abs(pc(~ind)).^2+abs(ps(~ind)).^2+abs(v(~ind)).^2),'b',...
            -sin(theta).*(abs(pc1).^2+abs(ps1).^2+abs(v1).^2), cos(theta).*(abs(pc1).^2+abs(ps1).^2+abs(v1).^2),'b',...
            sin(theta).*(abs(pc1).^2+abs(ps1).^2+abs(v1).^2), cos(theta).*(abs(pc1).^2+abs(ps1).^2+abs(v1).^2),'b',...
            -sin(theta(ind)),-cos(theta(ind)),'g',sin(theta(ind)),-cos(theta(ind)),'g',...
            [-12 12],[0 0],'k')
        axis image
        axis([-5 5 -2 2])
        text(3,1.5,['R.I = ' mnum2str(nm,1,3)])
        eval(['print -dpng -r200 RItmp' mint2str(j,2)])
    end
    
    % "c:\Program Files (x86)\ImageMagick-6.6.3-Q16\convert.exe"  -resize 50% RItmp*.png SiliconObjective.gif
end