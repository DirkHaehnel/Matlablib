function [para, facv, ind, taufit, wd] = AlexeyCavityFit(lifetimedata, spectrum, model, tau0, nn, flag)

if nargin<5 || isempty(nn)
    nn = 10;
end

if nargin<6 || isempty(flag)
    flag = [];
end

if isstr(lifetimedata)
    load(lifetimedata);
else
   lam = lifetimedata(:,1);
   tau = lifetimedata(:,2);
end

if isstr(spectrum)
    em = load(spectrum);
else
    em = spectrum;
end
em(:,2) = (em(:,2)-min(em(:,2)))/(max(em(:,2))-min(em(:,2)));

if isstr(model)
    load(model);
end

ind = 5:25;
wd = polyval(polyfit(model.tmax(ind),model.dv(ind),1),lam);
    
spec = abs(interp1(em(:,1),em(:,2),model.lambda,'linear',0));

fp = zeros(numel(model.dv),size(model.exc(1).fxc,2));
fv = fp;
for k=1:length(model.dv) 
    % model.excitation intensity
    fp(k,:) = model.exc(k).rho(:,1)'*(abs(model.exc(k).fxc(:,:,1)).^2+0.5*sum(abs(model.exc(k).fxc(:,:,2:end)).^2,3)+0.5*sum(abs(model.exc(k).fxs).^2,3) + ...
        abs(model.exc(k).fyc(:,:,1)).^2+0.5*sum(abs(model.exc(k).fyc(:,:,2:end)).^2,3)+0.5*sum(abs(model.exc(k).fys).^2,3))/2;
    fv(k,:) = model.exc(k).rho(:,1)'*(abs(model.exc(k).fzc(:,:,1)).^2+0.5*sum(abs(model.exc(k).fzc(:,:,2:end)).^2,3)+0.5*sum(abs(model.exc(k).fzs).^2,3));
end

% mpcolor(1:model.dz,dv,fp0.*detp0./sp0*asin(1.33/1.52)/200/(4/3*1.33)); colorbar
% mpcolor(1:model.dz,dv,fv0.*detv0./sv0*asin(1.33/1.52)/200/(4/3*1.33)); colorbar

facv = 1; %0.98:0.002:1.02;
taufit = zeros(numel(wd),numel(facv));
err = zeros(numel(facv),1);
c = err;
for jfac=1:numel(facv)    
    fac = facv(jfac);
    close;
    para(:,jfac) = Simplex('CavityFit',[0.8 1/1000],[0 0],[1 inf],[],[],fac*wd,tau,tau0,model.dv,fp,fv,model.bp,model.bv,model.lp,model.lv,spec,nn,flag);
    [err(jfac) taufit(:,jfac) c(jfac)] = CavityFit(para(:,jfac),fac*wd,tau,tau0,model.dv,fp,fv,model.bp,model.bv,model.lp,model.lv,spec,nn,flag);
end
ind = 1:numel(err);
ind = ind(err==min(err));
plot(lam,tau,'o',lam,taufit(:,ind))
s = {['\tau_0 = ' mnum2str(c(ind),1,2) ' ns'], ['\Phi_q = ' mnum2str(para(1,ind),1,2)], ['\tau_{rot} = ' mnum2str(1/6/para(2,ind),1,2)], ['\delta = ' mnum2str(1e2*facv(ind)-100,1,1) ' %']};
%s = {['\tau_0 = ' mnum2str(tau0,1,2) ' ns'], ['\Phi_q = ' mnum2str(para(1,ind),1,2)], ['\tau_{rot} = ' mnum2str(1/6/para(2,ind),1,2)]};
text(0.9*max(lam),0.9*max(tau),s);
ylabel('mean lifetime (ns)')
xlabel('max. transmission (nm)')

% phi = phi(ind);
% drot = drot(ind)
% delta = facv(ind);

