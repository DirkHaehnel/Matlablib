function [ex,ez,ex0,ez0] = DipoleField(zd,rhoa,za,n0,nd,ni,na,n1,d0,dd,di,da,d1)

% [ex,ez,ex0,ez0] = DipoleField(r,z,n0,n,n1,d0,d,d1) calculates the electric field amplitudes 
% at position r of a dipole at distance z from an interface within a layer 
% rhoa, za - position of field calculation
% zd  - molecule's distance from the bottom of its layer (na)
% n0 - vector of refractive indices of the stack below the molecule's layer
% na  - refracive index of the molecule's layer
% ni - refractive index of interveneing layers
% nd - refractive index of target layer
% n1 - vector of refractive indices of the stack above the the target's layer
% d* - layer thickness values

nmax = 1e3; kmax = 2e4; % most critical part: defines finesse of spacing in Simpson-rule integration
q = [(0.5:nmax)/nmax*nd nd+(0.5:kmax)/kmax*2e4*nd];
delta = [ones(1,nmax)*diff(q(1:2)) ones(1,kmax)*diff(q(end-1:end))];
wd = sqrt(nd^2 - q.^2);
wa = sqrt(na^2 - q.^2);

[~, ~, tpud, tsud] = Fresnel(wd,[nd ni na],di);
[rpud, rsud] = Fresnel(wd,[nd ni na n1],[di da d1]);
[rpdd, rsdd] = Fresnel(wd,[nd n0(end:-1:1)],d0(end:-1:1)); 
[rpua, rsua] = Fresnel(wa,[na n1],d1);
[rpda, rsda] = Fresnel(wa,[na ni(end:-1:1) nd n0(end:-1:1)],[di(end:-1:1) dd d0(end:-1:1)]);

tmp = exp(1i*wd*(dd-zd)+1i*wa*za);
tmp1 = (1 - rpdd.*exp(2*1i*wd*zd))./(1-rpud.*rpdd.*exp(2*1i*wd*dd))./(1-rpua.*rpda.*exp(2*1i*wa*da));
tmp1(isnan(tmp1)) = 1;
tp = tmp.*tpud.*tmp1;
tmp1 = (1 + rsdd.*exp(2*1i*wd*zd))./(1-rsud.*rsdd.*exp(2*1i*wd*dd))./(1-rsua.*rsda.*exp(2*1i*wa*da));
tmp1(isnan(tmp1)) = 1;
ts = tmp.*tsud.*tmp1;
tmp = exp(2*1i*wa*(da-za));
tmp(isnan(tmp)) = 0;
tp1 = tmp.*rpua.*tp;
ts1 = tmp.*rsua.*ts;

fac = 1i*q./wd.*delta;

b0 = besselj(0,q'*rhoa);
b1 = -1i*besselj(1,q'*rhoa);
b2 = -besselj(2,q'*rhoa);

% exact fields
dtotal = dd-zd + sum(di) + za;
rr = sqrt(dtotal^2+rhoa.^2);
ct = dtotal./rr; st = rhoa./rr;
ex0(1,:) = (nd^2./rr + nd*1i./rr.^2 - 1./rr.^3 - 0.5*(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*st.^2).*exp(1i*rr)/nd^2;
ex0(2,:) = (-0.5*(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*st.^2).*exp(1i*rr)/nd^2;
ex0(3,:) = (-0.5*(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*st.*ct).*exp(1i*rr)/nd^2;
ez0(1,:) = (nd^2./rr + nd*1i./rr.^2 - 1./rr.^3 -(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*ct.^2).*exp(1i*rr)/nd^2;
ez0(2,:) = (-(nd^2./rr+3*nd*1i./rr.^2-3./rr.^3).*st.*ct).*exp(1i*rr)/nd^2;

e0 = exp(1i*wd*(dd-zd+di+za));
tmp = fac.*(((tp-tp1).*wa/na-e0.*wd/nd).*wd/nd+(ts+ts1-e0));
%ex(1,:) = ex0(1,:) + 0.5*tmp(~isnan(tmp))*b0(~isnan(tmp),:); % x-component, cos(0*phi)
ex(1,:) = 0.5*tmp(~isnan(tmp))*b0(~isnan(tmp),:); % x-component, cos(0*phi)
tmp = fac.*(((tp-tp1).*wa/na-e0.*wd/nd).*wd/nd-(ts+ts1-e0));
%ex(2,:) = ex0(2,:) + 0.5*tmp(~isnan(tmp))*b2(~isnan(tmp),:);  % x-component, cos(2*phi) and y-component, sin(2*phi)
ex(2,:) = 0.5*tmp(~isnan(tmp))*b2(~isnan(tmp),:);  % x-component, cos(2*phi) and y-component, sin(2*phi)
tmp = fac.*((tp+tp1)/na-e0/nd).*q.*wd/nd;
%ex(3,:) = ex0(3,:) - 0.5*tmp(~isnan(tmp))*b1(~isnan(tmp),:); % z-component, cos(1*phi)
ex(3,:) = -0.5*tmp(~isnan(tmp))*b1(~isnan(tmp),:); % z-component, cos(1*phi)

tmp = exp(1i*wd*(dd-zd)+1i*wa*za);
tmp1 = (1 + rpdd.*exp(2*1i*wd*zd))./(1-rpud.*rpdd.*exp(2*1i*wd*dd))./(1-rpua.*rpda.*exp(2*1i*wa*da));
tmp1(isnan(tmp1)) = 1;
tp = tmp.*tpud.*tmp1;
tmp = exp(2*1i*wa*(da-za));
tmp(isnan(tmp)) = 0;
tp1 = tmp.*rpua.*tp;

tmp = fac.*((tp+tp1)/na-e0/nd).*q.^2/nd;
%ez(1,:) = ez0(1,:) + tmp(~isnan(tmp))*b0(~isnan(tmp),:); % z-component, cos(0*phi)
ez(1,:) = tmp(~isnan(tmp))*b0(~isnan(tmp),:); % z-component, cos(0*phi)
tmp = fac.*((tp-tp1).*wa/na-e0.*wd/nd).*q/nd;
%ez(2,:) = ez0(2,:) - tmp(~isnan(tmp))*b1(~isnan(tmp),:); % x-component, cos(1*phi) and y-component, sin(1*phi)
ez(2,:) = -tmp(~isnan(tmp))*b1(~isnan(tmp),:); % x-component, cos(1*phi) and y-component, sin(1*phi)


