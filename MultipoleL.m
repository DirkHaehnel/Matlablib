function [v,pc,ps,tp,ts,tp1,ts1,fac] = MultipoleL(theta,z,n1,n,n2,d1,d,d2,ord)

% [v,pc,ps] = MultipoleL(theta,z,n,d,ord) calculates the electric field amplitudes of multipole radiation 
% of order ord along emission angle theta 
% of a dipole at distance z from an interface within a layer 
% theta - direction of radiation downwards
% z  - molecule's distance from the bottom of its layer
% n1 - vector of refractive indices of the stack below the molecule's layer
% n  - refracive index of the molecule's layer
% n2 - vector of refractive indices of the stack above the the molecule's layer
% d1 - vector of layer thickness values of the stack below the molecule's layer ( length(d1)=length(n1)-1 )
% d  - thickness of molecule's layer
% d2 - vector of layer thickness values of the stack above the molecule's layer ( length(d2)=length(n2)-1 )
% multipole order

z = z(:)';
col = ones(size(z));
theta = abs(theta(:));
ind = theta<pi/2;

v = zeros(length(theta),length(z));
pc = v; ps = v;

if sum(ind)>0
    tmp = theta(ind);
    w = sqrt(n^2-n1(1)^2*sin(tmp).^2);
    [rpu, rsu, tpu, tsu] = Fresnel(w,[n n2],d2);
    [rpd, rsd, tpd, tsd] = Fresnel(w,[n n1(end:-1:1)],d1(end:-1:1));    
    tp = (tpd./(1-rpu.*rpd.*exp(2*i*w.'*d))).';
    ts = (tsd./(1-rsu.*rsd.*exp(2*i*w.'*d))).';
    tp1 = tp.*(rpu.*exp(2*i*w.'*d)).';
    ts1 = ts.*(rsu.*exp(2*i*w.'*d)).';
    fac = i^ord*sqrt(n1(1))*n1(1)*cos(tmp).*w.^(ord-1);
    ez = exp(i*w*z);
    v(ind,:) = ((fac.*n1(1)/n.*sin(tmp).*tp)*col).*ez + (-1)^ord*((fac.*n1(1)/n.*sin(tmp).*tp1)*col)./ez;
    pc(ind,:) = ((fac.*w/n.*tp)*col).*ez - (-1)^ord*((fac.*w/n.*tp1)*col)./ez;
    ps(ind,:) = ((fac.*ts)*col).*ez + (-1)^ord*((fac.*ts1)*col)./ez;
end
