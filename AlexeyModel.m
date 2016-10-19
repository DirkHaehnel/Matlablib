function model  = AlexeyModel(met, nwater, lamex, d0v, d1, dv)

if ischar(met)
    load metals
    ind = wavelength>=450;
    lambda = wavelength(ind);
else
    ind = met(:,1)>=450;
    lambda = met(ind,1);
end
    
nglass = 1.52;
NA = 1.2;
dz = 25;

if numel(d0v)==1
    d0v = d0v*ones(size(dv));
end


% emission calculation
lambda = lambda(1):5:lambda(end);
theta = (0.5:200)/200*asin(NA/nglass);
dtheta = mean(diff(theta));
bv = zeros(dz,numel(lambda),numel(dv));
bp = bv;
lv = bv;
lp = bv;
    
% emission water
n = nwater;
for jd = 1:length(dv)
    d = dv(jd);
    d0 = d0v(jd);
    for jw = 1:length(lambda)
        lamem = lambda(jw);
        if ischar(met)
            eval(['nagem = interp1(wavelength,' met ',lamem);']);
        else
            nagem = interp1(met(:,1),met(:,2),lamem);
        end
        n0 = [nglass nagem];
        n1 = [nagem nglass];
        zv = d/dz/2:d/dz:d;
        
        [~,~,~,~,qvd,qvu,qpd,qpu] = LifetimeL(zv/lamem*2*pi,n0,n,n1,d0/lamem*2*pi,d/lamem*2*pi,d1/lamem*2*pi);
        [v,pc,ps] = DipoleL(theta,zv/lamem*2*pi,n0,n,n1,d0/lamem*2*pi,d/lamem*2*pi,d1/lamem*2*pi);
        
        bv(:,jw,jd) = (sin(theta)*abs(v).^2)'*dtheta/(4/3*n);
        bp(:,jw,jd) = (0.5*sin(theta)*(abs(pc).^2+abs(ps).^2))'*dtheta/(4/3*n);
        lv(:,jw,jd) = 4/3*n./(qvd+qvu)';
        lp(:,jw,jd) = 4/3*n./(qpd+qpu)';
    end
    disp(jd)
end

% excitation water
if ischar(met)
    eval(['nagex = interp1(wavelength,' met ',lamex);']);
else
    nagex = interp1(met(:,1),met(:,2),lamex);
end
rhofield = [0 1];
fd = 3e3;
over = 4e3;
focpos = 0;
atf = [];
n0 = [nglass nagex];
n1 = [nagex nglass];
for jd = 1:length(dv)
    d = dv(jd);
    d0 = d0v(jd);
    resolution = [20 lamex/(d/dz)];
    exc(jd) = GaussExc(rhofield, [0 d/1e3], NA, fd, n0, n, n1, d0/1e3, d/1e3, d1/1e3, lamex/1e3, over, focpos, atf, resolution);
    disp(jd)
end

% j=11;FocusImage2D(model.exc(j).rho,model.exc(j).z,cat(4,cat(3,model.exc(j).fxc,model.exc(j).fxs),cat(3,model.exc(j).fyc,model.exc(j).fys)))
% j=11;FocusImage2D(model.exc(j).rho,model.exc(j).z,cat(3,model.exc(j).fzc,model.exc(j).fzs))

model.lambda = lambda;
%model.tmax = tmax;
model.dz = dz;
model.lv = lv;
model.lp = lp;
model.bv = bv;
model.bp = bp;
model.dv = dv;
model.d0v = d0v;
model.exc = exc;

% for j=30:5:40 model  = AlexeyModel('silver', 1.33, 488, j, 85, dv); eval(['save ModelSilver488ex' mint2str(j,2) 'nm']); end

% mpcolor(model.lambda,1e3*model.exc(1).z(1,:),model.lp(:,:,1)); colorbar
% mpcolor(model.lambda,1e3*model.exc(1).z(1,:),model.lv(:,:,1)); colorbar