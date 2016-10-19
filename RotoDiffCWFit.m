function [err, c, z] = RotoDiffCWFit(p, t, y, nflag, tt, alpha, ww)

% function RotoDiffCWFit(p, t, y, nflag, tt, alpha, ww) for fitting
% rotational diffusion FCS measurements
%
% input variables:
%
% p - fit parameter vector; its structure is determined by the nflag
% variable. If nflag is a simple number, than it is assumed that there are
% nflag independent rotors present. The number of iostropic rotors is then
% determined as n1 = 2*nflag(1) - length(p), and the number of symmetric top
% rotors as n2 = length(p) - sum(nflag). Correspondingly refer the first n1
% values of p to the radius values of the n1 isotropic rotors, and the next
% 2*n2 values of p to the radii values of the n2 symmetric top rotors. If
% nflag is a two-dimensional vector, then nflag(1) has the same meaning as
% befor, but additionally it is assumed that there are nflag(2)
% exponentials present in the signal, which are described by the
% p(n1+n2+1:end) exponential decay times.
%
% t - time axis
% y - four rotational diffusion FCS curves: xx-xx, xx-xy, xy-xx, xy-xy, 
% alpha - angle between excitation and emission dipole of dye (default = 0) 
% ww - beam waist radius of excitation focus (default = 0.4 mum)
% nflag - see description of p above (default = 1)
% tt - measurement temperature (default = 21°C)

t = t(:);
p = p(:)';

if nargin<6 || isempty(alpha)
    % alpha = interp1((ux(3,:)-uy(3,:))./(ux(3,:)+2*uy(3,:)),0:15:90,rr,'cubic');
    alpha = 0; % angle between excitation & emission dipole
end

if nargin<7 || isempty(ww)
    ww = inf; % focus radius in µm
end

if nargin<4 || isempty(nflag)
    nflag = [1 0];
end
if length(nflag)==1
    nflag = [nflag 0];
end

if nargin<5 || isempty(tt)
    tt = 21; % temperature
end

n1 = 2*nflag(1) - length(p) + nflag(2);
n2 = length(p) - sum(nflag);
zz = zeros(length(t),4*nflag(1));

if ww>0.611 || isinf(ww)
    alpha = alpha/180*pi;
    if n1>0
        for j=1:n1
            zz(:,4*(j-1)+1:4*j) = RotoDiffInf(1/Rad2RotoDiff(p(j),tt)/6,0,alpha,t);
        end
    end
    if n2>0
        for j=1:n2
            tau = 6*Rad2RotoDiff(p(n1+((j-1)*2+1:j*2)),tt);
            dd = 1./tau(1)-1./tau(2);
            dp = 1./tau(1);
            zz(:,4*n1+(4*(j-1)+1:4*j)) = RotoDiffInf(dp,dd,alpha,t);
        end
    end
else
    load ConfoAnisoCoefficients
    if ww<min(w0)
        ww = min(w0);
    end
    if n1>0
        for j=1:n1
            dp(j) = 1/Rad2RotoDiff(p(j),tt)/6;
        end
        for L=1:4
            tmpo = 0; tmop = 0; tmpp = 0; tmoo = 0;
            for M=0:L
                tmpo = tmpo + interp2(0:15:90, w0, squeeze(permute(po(L+1,M+1,:,:),[3 4 1 2])), alpha, ww, 'cubic');
                tmop = tmop + interp2(0:15:90, w0, squeeze(permute(op(L+1,M+1,:,:),[3 4 1 2])), alpha, ww, 'cubic');
                tmpp = tmpp + interp2(0:15:90, w0, squeeze(permute(pp(L+1,M+1,:,:),[3 4 1 2])), alpha, ww, 'cubic');
                tmoo = tmoo + interp2(0:15:90, w0, squeeze(permute(oo(L+1,M+1,:,:),[3 4 1 2])), alpha, ww, 'cubic');
            end
            for j=1:n1
                tmp = exp(-L*(L+1)*dp(j)*t);
                zz(:,4*(j-1)+1:4*j) = zz(:,4*(j-1)+1:4*j) + [tmpp*tmp, tmpo*tmp, tmop*tmp tmoo*tmp];
            end
        end
    end
    
    if n2>0
        for j=1:n2
            tau((j-1)*2+1:j*2) = 6*Rad2RotoDiff(p(n1+((j-1)*2+1:j*2)),tt);
        end
        dd = 1./tau(1:2:end)-1./tau(2:2:end);
        dp = 1./tau(1:2:end);
        for L=1:4
            for M=0:L
                tmpo = interp2(0:15:90, w0, squeeze(permute(po(L+1,M+1,:,:),[3 4 1 2])), alpha, ww, 'cubic');
                tmop = interp2(0:15:90, w0, squeeze(permute(op(L+1,M+1,:,:),[3 4 1 2])), alpha, ww, 'cubic');
                tmpp = interp2(0:15:90, w0, squeeze(permute(pp(L+1,M+1,:,:),[3 4 1 2])), alpha, ww, 'cubic');
                tmoo = interp2(0:15:90, w0, squeeze(permute(oo(L+1,M+1,:,:),[3 4 1 2])), alpha, ww, 'cubic');
                for j=1:n2
                    tmp = exp(-L*(L+1)*dp(j)*t-M^2*dd(j)*t);
                    zz(:,4*n1+(4*(j-1)+1:4*j)) = zz(:,4*n1+(4*(j-1)+1:4*j)) + [tmpp*tmp, tmpo*tmp, tmop*tmp tmoo*tmp];
                end
            end
        end
    end
end

ze = [];
if nflag(2)>0
    for j=n1+2*n2+1:length(p)
        %ze = [ze zz.*(exp(-t*(1./p(j)))*ones(1,size(zz,2)))];
        ze = [ze exp(-t*(1./p(j)))];
    end
end
    

if nargin>2 && ~isempty(y)
    for j=1:size(y,2)
        if nflag(2)==0
            c(:,j) = zz(:,j:4:end)\y(:,j);
            z(:,j) = zz(:,j:4:end)*c(:,j);
        else
            %c(:,j) = [zz(:,j:4:end) ze]\y(:,j);
            %z(:,j) = [zz(:,j:4:end) ze]*c(:,j);
            
            c(:,j) = lsqnonneg([zz(:,j:4:end) ze],y(:,j));
            z(:,j) = [zz(:,j:4:end) ze]*c(:,j);
%             p = Simplex('CommonFactorLSQ',0.1,0,[],[],[],[zz ze],y);
%             [~, z] = CommonFactorLSQ(p,[zz, ze],y);
        end
    end
    
    %plot(t,y*diag([1 6 3]),'o',t,z*diag([1 6 3])); drawnow
    subplot(4,1,1:3)
    %semilogx(t,y*diag([1 3 6]),'o',t,z*diag([1 3 6]));
    semilogx(t,y,'o',t,z);
    set(gca,'xtick',[]);
    ylabel('correlation (a.u.)')
    subplot(4,1,4)
    
    semilogx(t,(y-z));
    err = sum(sum((y-z).^2));
    ylabel('res.');
    xlabel('time (ns)');
    drawnow
else
    err = zz;
end

end

function zz = RotoDiffInf(dp,dd,gamma,t)
    zz(:,1) = (665 + 380*exp(12*dd*t) + 171*exp(16*dd*t) + 21600*exp(2*(6*dd + 7*dp)*t) + 7200*exp(2*(8*dd + 7*dp)*t) + 28*(35 + 20*exp(12*dd*t) + 9*exp(16*dd*t) + 1200*exp(2*(6*dd + 7*dp)*t) + 400*exp(2*(8*dd + 7*dp)*t))*cos(2*gamma) + (2835 + 1620*exp(12*dd*t) + 729*exp(16*dd*t) + 13920*exp(2*(6*dd + 7*dp)*t) + 4640*exp(2*(8*dd + 7*dp)*t))*cos(4*gamma))./(352800.*exp(4*(4*dd + 5*dp)*t));
    zz(:,2) = (-665 + 57*exp(16*dd*t) + 10080*exp(2*(6*dd + 7*dp)*t) + 4640*exp(2*(8*dd + 7*dp)*t) + 28*(-35 + 3*exp(16*dd*t) + 80*exp(2*(8*dd + 7*dp)*t))*cos(2*gamma) + (-2835 + 243*exp(16*dd*t) - 10080*exp(2*(6*dd + 7*dp)*t) + 800*exp(2*(8*dd + 7*dp)*t))*cos(4*gamma))./(352800.*exp(4*(4*dd + 5*dp)*t));
    zz(:,3) = (32*(35 + exp(16*dd*t) + 20*exp(2*(8*dd + 7*dp)*t))*cos(gamma)^4 - (1645 + 47*exp(16*dd*t) - 2560*exp(2*(8*dd + 7*dp)*t) - 420*exp(12*dd*t + 8*dp*t) + 5*(455 + 13*exp(16*dd*t) + 1568*exp(2*(6*dd + 7*dp)*t) + 288*exp(2*(8*dd + 7*dp)*t) - 364*exp(12*dd*t + 8*dp*t))*cos(2*gamma))*sin(gamma)^2)./(88200.*exp(4*(4*dd + 5*dp)*t));
    zz(:,4) = (665 - 380*exp(12*dd*t) + 171*exp(16*dd*t) - 21600*exp(2*(6*dd + 7*dp)*t) + 7200*exp(2*(8*dd + 7*dp)*t) + 28*(35 - 20*exp(12*dd*t) + 9*exp(16*dd*t) - 1200*exp(2*(6*dd + 7*dp)*t) + 400*exp(2*(8*dd + 7*dp)*t))*cos(2*gamma) + (2835 - 1620*exp(12*dd*t) + 729*exp(16*dd*t) - 13920*exp(2*(6*dd + 7*dp)*t) + 4640*exp(2*(8*dd + 7*dp)*t))*cos(4*gamma))./(352800.*exp(4*(4*dd + 5*dp)*t));
end

