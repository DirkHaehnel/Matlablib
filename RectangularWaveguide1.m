% program rectangular waveguide

clear all

lamem = 0.670; % wavelength in mum
width = 1.5/lamem*2*pi; % width of waveguide 
height = 0.5/lamem*2*pi; % height of waveguide 
nw = 1.51; % refractive index of waveguide
ng = 1.46; % refractive index of support
na = 1; % refractive index of surrounding medium

Kmax = 30;
wv = na+1e-3:1e-3:nw-1e-3; % axial wave vector values
%wv = 1.3:1e-3:1.35; % axial wave vector values

Kx = 20;
Ky = 20;
pos = [-width/2 + width*(0.5:Kx)/Kx, width/2*ones(1,Ky), width/2 - width*(0.5:Kx)/Kx, -width/2*ones(1,Ky); ...
    -height/2*ones(1,Kx), (-height/2+height*(0.5:Ky)/Ky), height/2*ones(1,Kx), (height/2-height*(0.5:Ky)/Ky)];

r = sqrt(sum(pos.^2));
phi = angle(pos(1,:)+i*pos(2,:));

r0 = sqrt((pos(1,:)-width/4).^2+(pos(2,:)-height/3).^2);
phi0 = angle((pos(1,:)-width/4)+i*(pos(2,:)-height/3));

tf = [-sin(phi(1:Kx)) cos(phi(Kx+1:Kx+Ky)) -sin(phi(Kx+Ky+1:2*Kx+Ky)) cos(phi(2*Kx+Ky+1:2*Kx+2*Ky))];
tr = [cos(phi(1:Kx)) sin(phi(Kx+1:Kx+Ky)) cos(phi(Kx+Ky+1:2*Kx+Ky)) sin(phi(2*Kx+Ky+1:2*Kx+2*Ky))];

clear Ein Eout Hin Hout
% coefficient vector: [s-in; p-in; s-out; p-out]
for jw = 1:length(wv)
    w = wv(jw);
    qin = sqrt(nw.^2-w.^2);
    qout = sqrt(na.^2-w.^2);    
    
    for m=-Kmax:Kmax
        ef = exp(i*m*phi);
        [fr, ff] = CylM(m,qin,r);
        Ein(Kmax+m+1,:) = [ef.*(ff.*tf+fr.*tr) zeros(size(ff))];
%        Hin(3*Kmax+m+2,:) = -i*nw^2*[ef.*(ff.*tf+fr.*tr) zeros(size(ff))];
%        [fr, ff] = CylM(m,qout,r,'h');
%        Eout(Kmax+m+1,:) = [ef.*(ff.*tf+fr.*tr) zeros(size(ff))];
%        Hout(3*Kmax+m+2,:) = -i*na^2*[ef.*(ff.*tf+fr.*tr) zeros(size(ff))];
        [fr, ff, fz] = CylN(m,qin,nw,r);
%        Ein(3*Kmax+m+2,:) = [ef.*(ff.*tf+fr.*tr) ef.*fz];
        Hin(Kmax+m+1,:) = -i*[ef.*(ff.*tf+fr.*tr) ef.*fz];
%        [fr, ff, fz] = CylN(m,qout,na,r,'h');
%        Eout(3*Kmax+m+2,:) = [ef.*(ff.*tf+fr.*tr) ef.*fz];
%        Hout(Kmax+m+1,:) = -i*[ef.*(ff.*tf+fr.*tr) ef.*fz];
    end 
    [fr, ff] = CylM(0,qin,r0);
    E0 = [ef.*(ff.*tf+fr.*tr) zeros(size(ff))];
    [fr, ff, fz] = CylN(0,qin,nw,r);
    H0 = -i*[ef.*(ff.*tf+fr.*tr) ef.*fz];
    
    %MM = [[Ein; -Eout] [Hin; -Hout]];
    MM = [Ein Hin];
%     weight = max(abs(MM),[],2);
%     MM = MM./(weight*ones(1,size(MM,2)));
    %tmp = real(cond(MM*MM.'))
    tmp = sum(abs([E0 H0]/MM))
    res(jw) = tmp(1);
end


