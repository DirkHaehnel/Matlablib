function [para, taufit, lam, tau, c] = AlexeyCavityFitFull(lifetimedata, spectrum, model, tau0, nn, flag)

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

spec = abs(interp1(em(:,1),em(:,2),model.lambda,'linear',0));

[~, ind] = intersect(model.tmax,lam);
model.tmax = model.tmax(ind);
model.lv = model.lv(:,:,ind);
model.lp = model.lp(:,:,ind);
model.bv = model.bv(:,:,ind);
model.bp = model.bp(:,:,ind);
model.dv = model.dv(ind);
model.d0v = model.d0v(ind);
model.exc = model.exc(ind);

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

close
para = Simplex('CavityFitFull',[0.8 tau0],[0 0],[1 inf],[],[],lam,tau,[],model.tmax,fp,fv,model.bp,model.bv,model.lp,model.lv,spec,nn,flag);
[~, taufit, c] = CavityFitFull(para,lam,tau,[],model.tmax,fp,fv,model.bp,model.bv,model.lp,model.lv,spec,nn,flag);
plot(lam,tau,'o',lam,taufit)
s = {['\tau_0 = ' mnum2str(c,1,2) ' ns'], ['\Phi_q = ' mnum2str(para(1),1,2)], ['\tau_{rot} = ' mnum2str(1/6/para(2),1,2) ' ns']};
text(0.9*max(lam),0.9*max(tau),s);
ylabel('mean lifetime (ns)')
xlabel('max. transmission (nm)')

% phi = phi(ind);
% drot = drot(ind)
% delta = facv(ind);

