function [v,pc,ps] = DipoleAniso(theta,z,n1,n,n2,d)

% [v,pc,ps] = DipoleAniso(theta,z,n1,n,n2,d) calculates the electric field amplitudes of a dipole 
% radiating out of an anisotropic (uniaxal) layer along emission angle theta 
% at a distance z away from the lower interface  
% theta - direction of radiation downwards
% z  - molecule's distance from the bottom of its layer
% n1 -  refractive index below the molecule's layer
% n  - birefringent refracive index of the molecule's layer (no,ne)
% n2 - refractive index above the the molecule's layer
% d  - thickness of molecule's layer


z = z(:)';
col = ones(size(z));
theta = abs(theta(:));
ind = theta<pi/2;

v = zeros(length(theta),length(z));
pc = v; ps = v;

no = n(1);
ne = n(2);
if sum(ind)>0
    tmp = theta(ind);
    q = n1*sin(tmp);
    wo = sqrt(no^2-q.^2);
    we = no*sqrt(1-q.^2/ne^2);    
    [rpu, rsu, tpu, tsu] = FresnelAniso(q, n, [n2 n2]);
    [rpd, rsd, tpd, tsd] = FresnelAniso(q, n, [n1 n1]);    
    tp = tpd./(1-rpu.*rpd.*exp(2*i*we*d));
    ts = tsd./(1-rsu.*rsd.*exp(2*i*wo*d));
    tp1 = tp.*rpu.*exp(2*i*we*d);
    ts1 = ts.*rsu.*exp(2*i*wo*d);
    rt = sqrt(we.^2/no^4+q.^2/ne^4);
    fac = sqrt(n1)*n1^2*cos(tmp);
    ezo = exp(i*wo*z);
    eze = exp(i*we*z);
    ezo1 = conj(ezo);
    eze1 = conj(eze);    
    v(ind,:) = ((fac.*q.*no^2/ne^2./we.*rt.*tp)*col).*eze + ...
        ((fac.*q.*no^2/ne^2./we.*rt.*tp1)*col).*eze1;
    pc(ind,:) = ((fac.*rt.*tp)*col).*eze - ((fac.*rt.*tp1)*col).*eze1;
    ps(ind,:) = ((fac./wo.*ts)*col).*ezo + ((fac./wo.*ts1)*col).*ezo1;
end
