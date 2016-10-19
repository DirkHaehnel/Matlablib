% program rectangular waveguide

close

lamem = 0.670; % wavelength in mum
width = pi; % 1.00/lamem*2*pi; % width of waveguide 
height = pi/2; %0.2500/lamem*2*pi; % height of waveguide 
nw = 1.51; % refractive index of waveguide
ng = 1.46; % refractive index of support
na = 1; % refractive index of surrounding medium

Kmax = 1e4;
psi = pi/2*(0.5:Kmax)'/Kmax; % azimuthal angles of considered plane wave modes
qv = na + (1:1e3-1)/1e3*(ng-na); % axial wave vector values

% coeffieicnt vector: [s-in; p-in; s-out; p-out]
for jq = 1:length(qv)
    q = qv(jq);
    kvec = [cos(psi) sin(psi)]*sqrt(nw^2-q^2); % wave vector
    rx = Fresnel(kvec(:,1),nw,na);
    ry = Fresnel(kvec(:,2),nw,na);
    plot(psi,mod(angle(rx.^2)+2*kvec(:,1)'*width,2*pi)+mod(angle(ry.^2)+2*kvec(:,2)'*height,2*pi)); 
    axis([0 pi/2 0 20])
    drawnow
end


