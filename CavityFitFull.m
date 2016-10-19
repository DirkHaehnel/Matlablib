function [err,tau,c] = CavityFitFull(p,lam,taud,tau0,tmax,fp,fv,bp,bv,lp,lv,spec,nn,flag)

taud = taud(:);
if nargin<13 || isempty(nn)
    nn = 10;
end
a0 = 0.5;
mm = 2*(0:nn-1)';

for k=1:length(lam)
    ind = tmax==lam(k);
    
    % rel. detection efficiency 
    dpp(k,:) = spec*(bp(:,:,ind).*((1-p(1))*lp(:,:,ind)+p(1)))'/sum(spec);
    dvv(k,:) = spec*(bv(:,:,ind).*((1-p(1))*lv(:,:,ind)+p(1)))'/sum(spec);
    dpv(k,:) = spec*(bp(:,:,ind).*((1-p(1))*lp(:,:,ind)+p(1).*lp(:,:,ind)./lv(:,:,ind)) + bv(:,:,ind).*((1-p(1))*lv(:,:,ind)+p(1)*lv(:,:,ind)./lp(:,:,ind)))'/sum(spec);
    
    % emission rate
    sp(k,:) = spec*(1./lp(:,:,ind)')/sum(spec);
    sv(k,:) = spec*(1./lv(:,:,ind)')/sum(spec);
end
ds = p(1)*(sv-sp); %/tau0;

tau = zeros(size(fp,1),1);
for j=1:size(fp,1)
    nrm = 0;
    for k=1:size(fp,2)
        % initial Legendre coefficent
        a2 = (fv(j,k)-fp(j,k))/(fv(j,k)+2*fp(j,k));
        
        % rotational diffusion matrix
        %[V,D] = eig(diag(ds(j,k)*mm(2:end).*(mm(2:end)-1)./(2*mm(2:end)-3)./(2*mm(2:end)-1),-1) + diag((1-p(1)+p(1)*sp(j,k))/tau0 + ds(j,k)*(2*mm.*(mm+1)-1)./(2*mm-1)./(2*mm+3) + mm.*(mm+1)*p(2)) + diag(ds(j,k)*(mm(1:end-1)+1).*(mm(1:end-1)+2)./(2*mm(1:end-1)+3)./(2*mm(1:end-1)+5),1));
        [V,D] = eig(diag(ds(j,k)*mm(2:end).*(mm(2:end)-1)./(2*mm(2:end)-3)./(2*mm(2:end)-1),-1) + diag((1-p(1)+p(1)*sp(j,k)) + ds(j,k)*(2*mm.*(mm+1)-1)./(2*mm-1)./(2*mm+3) + mm.*(mm+1)*p(2)) + diag(ds(j,k)*(mm(1:end-1)+1).*(mm(1:end-1)+2)./(2*mm(1:end-1)+3)./(2*mm(1:end-1)+5),1));
        Vi = inv(V);
        
        for r=1:nn
            ar = a0*Vi(r,1) + a2*Vi(r,2); % projection of initial state on eigenvectors
            tau(j) = tau(j) + ar*(2/15*(3*dvv(j,k)+2*dpv(j,k)+8*dpp(j,k))*V(1,r) + ...
                4/105*(6*dvv(j,k)+dpv(j,k)-8*dpp(j,k))*V(2,r) + ...
                16/315*(dvv(j,k)-dpv(j,k)+dpp(j,k))*V(3,r))/D(r,r)^2; % integral <t> F(t)
            nrm = nrm + ar*(2/15*(3*dvv(j,k)+2*dpv(j,k)+8*dpp(j,k))*V(1,r) + ...
                4/105*(6*dvv(j,k)+dpv(j,k)-8*dpp(j,k))*V(2,r) + ...
                16/315*(dvv(j,k)-dpv(j,k)+dpp(j,k))*V(3,r))/D(r,r); % integral F(t)
        end
    end
    tau(j) = tau(j)/nrm;
end

if nargin>13 && ~isempty(flag)
    paffine = Simplex('AffineFit',[0 1],[],[],[],[],wd,taud,dv,tau,[],[],-1,1)
    [err, c, tau] = AffineFit(paffine,wd,taud,dv,tau,[],[],-1,1);
else
    if isempty(tau0)
        c = lsqnonneg(tau,taud);
        tau = c*tau;
    else
        tau = tau0*tau;
        c = tau0;
    end
end

plot(lam,taud,'o',lam,tau); drawnow
err = sum(abs((taud-tau).^2));

