% figure for astigmatism paper

flag = 1;

close all
NA = 1.2;
chimax = asin(NA/1.333);
chi = (0:1e2)/100*chimax;
phi = pi/2 + (0:100)/200*2*pi;
foc = 200;
len = 2.1*foc*sin(chimax);

switch flag
    case 1 % perfect image 
        d = 1.5;
        surf([0 len], [-len len], -foc*ones(2,2), zeros(2,2),'facealpha', 0.1) 
        hold on
        surf([0 0; len len], -len*ones(2,2), -foc*[1 d; 1 d], zeros(2,2),'facealpha', 0.3) 
        surf(zeros(2,2), [-len -len; len len], -foc*[1 d; 1 d], zeros(2,2),'facealpha', 0.2)
        plot3([0 0], [0 0], [-d*foc 0.5*foc], 'linewidth', 1, 'color', [0,0,0])
        theta0 = 0.8*chimax;
        rad = foc*tan(theta0);
        theta = asin(1.333*sin(theta0)/1.51);
        surf([0; -rad]*cos(phi),[0; -rad]*sin(phi),[0; -foc]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', 0.4)
        surf([-rad; -rad-tan(theta)*(d-1)*foc]*cos(phi),[-rad; -rad-tan(theta)*(d-1)*foc]*sin(phi),[-foc; -d*foc]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', 0.4)
        len = 120;
        surf(-len*sin(chi')*cos(phi),len*sin(chi')*sin(phi),-len*cos(chi')*ones(size(phi)),0.*ones(length(chi), length(phi)), 'facealpha', 0.1)
        plot3([0 0],[-rad, -rad],-foc+[-100 100], 'linewidth', 1, 'color', [0 0 0]);
        plot3(0*chi(chi<theta0),-rad+90*sin(chi(chi<theta0)),-foc+90*cos(chi(chi<theta0)),'linewidth', 1, 'color', [0 0 0])
        plot3(0*chi(chi<theta),-rad-90*sin(chi(chi<theta)),-foc-90*cos(chi(chi<theta)),'linewidth', 1, 'color', [0 0 0])
        [x,y,z] = sphere(50);
        surf(3*x,3*y,3*z,zeros(size(x)));
        hold off
        
    case 2 % coverslide 
        d = 1.5;
        dd = 100;
        surf([0 len], [-len len], dd-foc*ones(2,2), zeros(2,2),'facealpha', 0.1) 
        hold on
        surf([0 len], [-len len], -foc*ones(2,2), zeros(2,2),'facealpha', 0.1)         
        surf([0 0; len len], -len*ones(2,2), [dd 0; dd 0]-foc*[1 d; 1 d], zeros(2,2),'facealpha', 0.3) 
        surf(zeros(2,2), [-len -len; len len], [dd 0; dd 0]-foc*[1 d; 1 d], zeros(2,2),'facealpha', 0.2)
        plot3([0 0], [0 0], [-d*foc 0.5*foc], 'linewidth', 1, 'color', [0,0,0])
        theta0 = 0.8*chimax;
        theta = asin(1.333*sin(theta0)/1.51);
        surf([0; -tan(theta0)*foc+tan(theta)*dd]*cos(phi),...
            [0; -tan(theta0)*foc+tan(theta)*dd]*sin(phi),...
            [dd-tan(theta)/tan(theta0)*dd; -foc+dd]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', 0.4)
        surf([-tan(theta0)*foc+tan(theta)*dd; -tan(theta0)*foc-tan(theta)*(d-1)*foc]*cos(phi),...
            [-tan(theta0)*foc+tan(theta)*dd; -tan(theta0)*foc-tan(theta)*(d-1)*foc]*sin(phi),...
            [dd-foc; -d*foc]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', 0.4)
        theta0 = 0.4*chimax;
        theta = asin(1.333*sin(theta0)/1.51);
        surf([0; -tan(theta0)*foc+tan(theta)*dd]*cos(phi),...
            [0; -tan(theta0)*foc+tan(theta)*dd]*sin(phi),...
            [dd-tan(theta)/tan(theta0)*dd; -foc+dd]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', 0.4)
        surf([-tan(theta0)*foc+tan(theta)*dd; -tan(theta0)*foc-tan(theta)*(d-1)*foc]*cos(phi),...
            [-tan(theta0)*foc+tan(theta)*dd; -tan(theta0)*foc-tan(theta)*(d-1)*foc]*sin(phi),...
            [dd-foc; -d*foc]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', 0.4)
        
        len = 120;
        %    surf(-len*sin(chi')*cos(phi),len*sin(chi')*sin(phi),-len*cos(chi')*ones(size(phi)),0.*ones(length(chi), length(phi)), 'facealpha', 0.1)
        [x,y,z] = sphere(50);
        surf(3*x,3*y,3*z,zeros(size(x)));
        hold off
        
    case 3 % astigmatism 
        d = 1.5;
        dd = 25;
        surf([0 len], [-len len], -foc*ones(2,2), zeros(2,2),'facealpha', 0.1) 
        hold on
        surf([0 0; len len], -len*ones(2,2), -foc*[1 d; 1 d], zeros(2,2),'facealpha', 0.3) 
        surf(zeros(2,2), [-len -len; len len], -foc*[1 d; 1 d], zeros(2,2),'facealpha', 0.2)
        plot3([0 0], [0 0], [-d*foc 0.5*foc], 'linewidth', 1, 'color', [0,0,0])
        theta0 = 0.8*chimax;
        theta = asin(1.333*sin(theta0)/1.51);
        surf([zeros(1,length(phi)); tan(theta0)*(-foc+dd*cos(2*phi)).*cos(phi)],...
            [zeros(1,length(phi)); tan(theta0)*(-foc+dd*cos(2*phi)).*sin(phi)],...
            [dd*cos(2*phi); -foc*ones(1,length(phi))],0.6*ones(2,length(phi)), 'facealpha', 0.4)
        surf([tan(theta0)*(-foc+dd*cos(2*phi)).*cos(phi); (tan(theta0)*(-foc+dd*cos(2*phi))-tan(theta)*(d-1)*foc).*cos(phi)],...
            [tan(theta0)*(-foc+dd*cos(2*phi)).*sin(phi); (tan(theta0)*(-foc+dd*cos(2*phi))-tan(theta)*(d-1)*foc).*sin(phi)],...
            [-foc; -d*foc]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', 0.4)
        len = 120;
        surf(-len*sin(chi')*cos(phi),len*sin(chi')*sin(phi),-len*cos(chi')*ones(size(phi)),0.*ones(length(chi), length(phi)), 'facealpha', 0.1)
        surf(-len*sin(chi')*cos(phi),len*sin(chi')*sin(phi),-sin(chi').^2*dd*cos(2*phi)-len*cos(chi')*ones(size(phi)),0.*ones(length(chi), length(phi)), 'facealpha', 0.3)
        [x,y,z] = sphere(50);
        surf(3*x,3*y,3*z,zeros(size(x)));
        hold off
        
    case 4 % refractive index 
        d = 1.5;
        surf([0 len], [-len len], -foc*ones(2,2), zeros(2,2),'facealpha', 0.1) 
        hold on
        surf([0 0; len len], -len*ones(2,2), -foc*[1 d; 1 d], zeros(2,2),'facealpha', 0.3) 
        surf(zeros(2,2), [-len -len; len len], -foc*[1 d; 1 d], zeros(2,2),'facealpha', 0.2)
        plot3([0 0], [0 0], [-d*foc 0.5*foc], 'linewidth', 1, 'color', [0,0,0])
        theta0 = 0.8*chimax;
        rad = tan(theta0)*foc;
        theta = asin(1.333*sin(theta0)/1.51);
        theta0 = asin(1.51*sin(theta)/1.4);
        surf([-rad; -rad-tan(theta)*(d-1)*foc]*cos(phi),...
            [-rad; -rad-tan(theta)*(d-1)*foc]*sin(phi),...
            [-foc; -d*foc]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', 0.4)
        surf([0; -rad]*cos(phi), [0; -rad]*sin(phi),...
            [-foc+rad/tan(theta0); -foc]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', 0.4)
        theta0 = 0.4*chimax;
        rad = tan(theta0)*foc;
        theta = asin(1.333*sin(theta0)/1.51);
        theta0 = asin(1.51*sin(theta)/1.4);
        surf([-rad; -rad-tan(theta)*(d-1)*foc]*cos(phi),...
            [-rad; -rad-tan(theta)*(d-1)*foc]*sin(phi),...
            [-foc; -d*foc]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', 0.4)
        surf([0; -rad]*cos(phi), [0; -rad]*sin(phi),...
            [-foc+rad/tan(theta0); -foc]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', 0.4)
        [x,y,z] = sphere(50);
        surf(3*x,3*y,3*z,zeros(size(x)));
        hold off
end


caxis([0 1])
shading interp
axis image
axis off
view([-58 11])
