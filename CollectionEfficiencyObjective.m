% program CollectionEfficiencyObjective for calculating coll. eff. of an objective
% with num. apert. NA

nm = 1.333; % sample ref. index
ng = 1.523; % ref. index immersion medium of objective
NA = ng;%ng*sin(61.09/180*pi); % num. ap.
zv = 0:pi/50:5*pi; % vector of considered distances between molecules and surface

% calculate max. collection angle
theta_max = asin(NA/ng);

dtheta = pi/1e3;
theta = dtheta/2:dtheta:pi;
ind = theta<=theta_max;

for j=1:length(zv)

    % calculate angular dist. of emission into glass
    [v, pc, ps] = DipoleL(theta,zv(j),ng,nm,nm,[],zv(j),[]);
    
    % calculate angular dist. of emission into sample
    [v1,pc1,ps1] = DipoleL(theta,0,nm,nm,ng,[],zv(j),[]);
  
    % calculate collection efficiency
    cefv(j) = sum(sin(theta(ind))*abs(v(ind)).^2)/(sum(sin(theta)*abs(v).^2)+sum(sin(theta)*abs(v1).^2));
    cefp(j) = sum(sin(theta(ind))*(abs(pc(ind)).^2+abs(ps(ind)).^2))/...
        (sum(sin(theta)*(abs(pc).^2+abs(ps).^2))+sum(sin(theta)*(abs(pc1).^2+abs(ps1).^2)));    
    
end

plot(zv/pi/2,cefv,zv/pi/2,cefp)
xlabel('distance (\lambda)');
ylabel('collection efficiency');
legend({'vert. dipole','hor. dipole'})

