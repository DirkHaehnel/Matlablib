function [ex, ez] = DipoleFree(xx,yy,zz,x,y,z,nd)

if nargin<6 || isempty(nd)
    nd = 1;
end

rr = sqrt((xx-x).^2+(yy-y).^2+(zz-z).^2);
theta = angle((zz-z) + 1i*sqrt((xx-x).^2+(yy-y).^2));
ct = cos(theta); st = sin(theta);

ex(1,:,:,:) = (nd^2./rr + nd*1i./rr.^2 - 1./rr.^3 - 0.5*(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*st.^2).*exp(1i*rr)/nd^2;  % x-component, cos(0*phi)
ex(2,:,:,:) = (-0.5*(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*st.^2).*exp(1i*rr)/nd^2; % x-component, cos(2*phi) and y-component, sin(2*phi)
ex(3,:,:,:) = (-0.5*(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*st.*ct).*exp(1i*rr)/nd^2; % z-component, cos(1*phi)
ez(1,:,:,:) = (nd^2./rr + nd*1i./rr.^2 - 1./rr.^3 -(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*ct.^2).*exp(1i*rr)/nd^2; % z-component, cos(0*phi)
ez(2,:,:,:) = (-(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*st.*ct).*exp(1i*rr)/nd^2; % x-component, cos(1*phi) and y-component, sin(1*phi)


