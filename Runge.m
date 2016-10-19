if 1
    max_psi = 5e2;
    max_sig = 3e2;
    
    gpsi = zeros(max_psi,max_sig);
    psi = (0.5:max_psi)'/max_psi*2*pi;
    sig = (0:max_sig-1)/(max_sig-1);
    
	dsig = mean(diff(sig));
    dpsi = mean(diff(psi));
    
    gpsi = ones(max_psi,1)*((0:max_sig-1)/999*(-log(0.1)) + log(0.1));
    gpsi(:,end) = log(1 + 5*cos(3*psi).^2);
    
    nn = max_psi*max_sig;
    e = ones(nn,1);
%     M = spdiags([e e -4*e e e],[-max_psi -1:1 max_psi],nn,nn);
%     for j=1:max_sig
%         M((j-1)*max_psi+1, j*max_psi) = 1;
%         M(j*max_psi, (j-1)*max_psi+1) = 1;
%         if j<max_sig
%             M(j*max_psi, j*max_psi+1) = 0;
%             M(j*max_psi+1, j*max_psi) = 0;
%         end
%     end

%     tmp = - M(max_psi+1:end-max_psi,max_psi+1:end-max_psi)\(M(max_psi+1:end-max_psi,1:max_psi)*gpsi(:,1) + M(max_psi+1:end-max_psi,end-max_psi+1:end)*gpsi(:,end));
%     gpsi(:,2:end-1) = reshape(tmp,max_psi,max_sig-2);
%     gpsi = exp(gpsi);

    Mpsi = spdiags([e -2*e e],[-1:1],nn,nn);
    Msig = spdiags([e -2*e e],[-max_psi 0 max_psi],nn,nn);
    for j=1:max_sig
        Msig((j-1)*max_psi+1, j*max_psi) = 1;
        Msig(j*max_psi, (j-1)*max_psi+1) = 1;
        if j<max_sig
            Mpsi(j*max_psi, j*max_psi+1) = 0;
            Mpsi(j*max_psi+1, j*max_psi) = 0;
        end
    end

    tmp = - (Msig(max_psi+1:end-max_psi,max_psi+1:end-max_psi)+Mpsi(max_psi+1:end-max_psi,max_psi+1:end-max_psi))\...
        ((Msig(max_psi+1:end-max_psi,1:max_psi)+Mpsi(max_psi+1:end-max_psi,1:max_psi))*gpsi(:,1) + ...
        (Msig(max_psi+1:end-max_psi,end-max_psi+1:end)+Mpsi(max_psi+1:end-max_psi,end-max_psi+1:end))*gpsi(:,end));
    gpsi(:,2:end-1) = reshape(tmp,max_psi,max_sig-2);
    gpsi = exp(gpsi);

    nshift = zeros(max_psi,1);
    ffun = nshift;
end

if 1
    close all
    tmax = 1e3;
    dt = 1e-3;
    resx = zeros(length(psi),tmax);
    resy = resx;
    for t = 1:tmax
        
        [fsig,fpsi] = gradient(-1./gpsi,dsig,dpsi);
        fpsi = fpsi*dsig.*gpsi.^2;
        p = zeros(length(psi),length(sig),2);
        for j=1:max_psi
            p(j,1,1) = 0.1*cos(psi(j));
            p(j,1,2) = 0.1*sin(psi(j));
            v = [cos(psi(j)) sin(psi(j))];
            p(j,2,1) = p(j,1,1) + v(1)*dsig*gpsi(j,1);
            p(j,2,2) = p(j,1,2) + v(2)*dsig*gpsi(j,1);
            for k=2:length(sig)-1
                v = sqrt(1-fpsi(j,k)^2)*v - fpsi(j,k)*[-v(2) v(1)];
                v = v/sqrt(sum(v.^2));
                p(j,k+1,1) = p(j,k,1) + v(1)*dsig*gpsi(j,k);
                p(j,k+1,2) = p(j,k,2) + v(2)*dsig*gpsi(j,k);
            end
        end
        plot(squeeze(p([1:5:end 1],:,1))',squeeze(p([1:5:end 1],:,2))','b',squeeze(p([1:end 1],[1:5:end end],1)),squeeze(p([1:end 1],[1:5:end end],2)),'r');
        axis image
        title(t)
        drawnow
        
        resx(:,t) = squeeze(p(:,end,1));
        resy(:,t) = squeeze(p(:,end,2));
        
        npsi = gradient(nshift,dpsi);
        v = 1./(sum(gpsi,2)*dsig);
        vpsi = gradient(v.^2,dpsi);
        nshift = nshift - dt*vpsi./gpsi(:,end)/4;
        gpsi(:,end) = gpsi(:,end) + dt*(v.*fsig(:,end) + npsi./gpsi(:,end));
        plot(1:max_psi,v.^2,1:max_psi,nshift./gpsi(:,end)*100,1:max_psi,v.*fsig(:,end)/10)
        
        
        gpsi = log(gpsi);
        tmp = - M(max_psi+1:end-max_psi,max_psi+1:end-max_psi)\(M(max_psi+1:end-max_psi,1:max_psi)*gpsi(:,1) + M(max_psi+1:end-max_psi,end-max_psi+1:end)*gpsi(:,end));
        gpsi(:,2:end-1) = reshape(tmp,max_psi,max_sig-2);
        gpsi = exp(gpsi);

    end
end

return 

plot(squeeze(p([1:5:end/5 end/5],:,1))',squeeze(p([1:5:end/5 end/5],:,2))','b',squeeze(p(1:end/5,[1:5:end end],1)),squeeze(p(1:end/5,[1:5:end end],2)),'r');
axis image
axis off

plot(squeeze(p([1:5:end/10 end/10],:,1))',squeeze(p([1:5:end/10 end/10],:,2))','b',squeeze(p(1:end/10,[1:5:end end],1)),squeeze(p(1:end/10,[1:5:end end],2)),'r');
axis image
axis off