function [err, c, zz, z] = CumulantFitFun(p, t, y, pic)

% p - single value decay time or vector of decay times
% t - time axis
% y - decay curve
% pic - whether to draw figures during fitting
% pic = 1 - linear plot with residuals
% pic = 2 - semilogy plot
% pic > 2 - semilogx plot
% nrm - whether to normalize fit functions or not (default not)
% pos - if empty, fitting is done with lsq, otherwise with lsqnonneg

if length(t)<size(y,1)
    c = t;
    t = y(:);
    t = t(:);
    m = size(y,1);
    zz = [ones(m,1) exp(-abs(t)*p) exp(-2*abs(t)*p) exp(-3*abs(t)*p) exp(-4*abs(t)*p)];
    err(:,1) = zz(:,[1 2])*c(:,1);
    err(:,2) = zz(:,[1 3])*c(:,2);
    err(:,3) = zz(:,[1 4 5])*c(:,3);
else
    t = t(:);
    m = size(y,1);
    zz = [ones(m,1) exp(-abs(t)*p) exp(-2*abs(t)*p) exp(-3*abs(t)*p) exp(-4*abs(t)*p)];

    c(:,1) = lsqnonneg(zz(:,[1 2]),y(:,1));
    z(:,1) = zz(:,[1 2])*c(:,1);
    c(:,2) = lsqnonneg(zz(:,[1 3]),y(:,2));
    z(:,2) = zz(:,[1 3])*c(:,2);
    c(:,3) = lsqnonneg(zz(:,[1 4 5]),y(:,3));
    z(:,3) = zz(:,[1 4 5])*c(:,3);
   
    if nargin>3 && ~isempty(pic)
        if pic==1
            subplot(4,1,1:3)
            plot(tt, y, 'ob', tt, z, 'r'); 
            axis tight
            subplot(4,1,4)
            plot(tt, (y-z)./sqrt(abs(z)),t,0*t)
            axis tight
            drawnow
        elseif pic==2
            semilogy(tt, y, 'o', tt, z); drawnow
        else
            semilogx(tt, y, 'o', tt, z); drawnow
        end
    end

    %err = sum(sum((y-z).^2./abs(z)));
    %err = sum((y-z).^2);
    ind = y>0;
    err = sum(y(ind).*log(y(ind)./z(ind))-y(ind)+z(ind));
end

