% program rectangular waveguide

lamem = 0.670; % wavelength in mum
width = 3.500/lamem*2*pi; % width of waveguide 
height = 3.500/lamem*2*pi; % height of waveguide 
nw = 1.51; % refractive index of waveguide
ng = 1.46; % refractive index of support
na = 1; % refractive index of surrounding medium

Kmax = 80;
psi = 2*pi*(0.5:Kmax)'/Kmax; % azimuthal angles of considered plane wave modes
qv = (0:1e3-1)/1e3*nw; % axial wave vector values

Kx = Kmax/4;
Ky = Kmax/4;
pos = [-width/2 + width*(0.5:Kx)/Kx, width/2*ones(1,Ky), width/2 - width*(0.5:Kx)/Kx, -width/2*ones(1,Ky); ...
    -height/2*ones(1,Kx), (-height/2+height*(0.5:Ky)/Ky), height/2*ones(1,Kx), (height/2-height*(0.5:Ky)/Ky)];

% coeffieicnt vector: [s-in; p-in; s-out; p-out]
for jq = 1:length(qv)
    q = qv(jq);
    kvec = [[cos(psi) sin(psi)]*sqrt(nw^2-q^2) q*ones(size(psi))]; % wave vector
    ps = [-kvec(:,2) kvec(:,1) zeros(size(psi))]; ps = ps./sqrt(sum(ps.^2,2)*[1 1 1]);
    pt = cross(kvec,ps); pt = pt./sqrt(sum(pt.^2,2)*[1 1 1]);
    
    fac = exp(i*(kvec(:,1)*pos(1,:) + kvec(:,2)*pos(2,:)));
    Ein = [[ps(:,1)*ones(1,Kx) ps(:,2)*ones(1,Ky) ps(:,1)*ones(1,Kx) ps(:,2)*ones(1,Ky)].*fac (ps(:,3)*ones(1,size(pos,2))).*fac];
    Ein = [Ein; [[pt(:,1)*ones(1,Kx) pt(:,2)*ones(1,Ky) pt(:,1)*ones(1,Kx) pt(:,2)*ones(1,Ky)].*fac (pt(:,3)*ones(1,size(pos,2))).*fac]];    
    Hin = [[pt(:,1)*ones(1,Kx) pt(:,2)*ones(1,Ky) pt(:,1)*ones(1,Kx) pt(:,2)*ones(1,Ky)].*fac (pt(:,3)*ones(1,size(pos,2))).*fac];
    Hin = nw*[Hin; -[[ps(:,1)*ones(1,Kx) ps(:,2)*ones(1,Ky) ps(:,1)*ones(1,Kx) ps(:,2)*ones(1,Ky)].*fac (ps(:,3)*ones(1,size(pos,2))).*fac]];    

	kvec = [[cos(psi) sin(psi)]*sqrt(na^2-q^2)]; % wave vector
    fac = exp(i*(kvec(:,1)*pos(1,:) + kvec(:,2)*pos(2,:)));
    fac(abs(fac)>1) = 0;
    Eout = [[ps(:,1)*ones(1,Kx) ps(:,2)*ones(1,Ky) ps(:,1)*ones(1,Kx) ps(:,2)*ones(1,Ky)].*fac (ps(:,3)*ones(1,size(pos,2))).*fac];
    Eout = [Eout; [[pt(:,1)*ones(1,Kx) pt(:,2)*ones(1,Ky) pt(:,1)*ones(1,Kx) pt(:,2)*ones(1,Ky)].*fac (pt(:,3)*ones(1,size(pos,2))).*fac]];    
    Hout = [[pt(:,1)*ones(1,Kx) pt(:,2)*ones(1,Ky) pt(:,1)*ones(1,Kx) pt(:,2)*ones(1,Ky)].*fac (pt(:,3)*ones(1,size(pos,2))).*fac];
    Hout = na*[Hout; -[[ps(:,1)*ones(1,Kx) ps(:,2)*ones(1,Ky) ps(:,1)*ones(1,Kx) ps(:,2)*ones(1,Ky)].*fac (ps(:,3)*ones(1,size(pos,2))).*fac]];    
    
    MM = [[Ein; -Eout] [Hin; -Hout]];
    res(jq) = det(MM.');
end


