function [err, c, z] = Rigler_Cos(p, t, y, m, bld)

% Function [err, c, z] = Rigler_Cos(p, t, y) for FCS data with harmonic modulation. 
% Calculates the error err, amplitudes c, and fit curve z for the 
% FCS data y with the model offset + 1/(p(1)+t)/sqrt(p(1)/p(2)^2+t)*(1 + cos(w t) + exp(-t/p(3)) + ...)
% p(1) : diffusion time
% p(2) : aspect ratio a
% p(3:end) : m modulation terms and the rest inverse of tripplet lifetimes
% m : number of modulation terms

    t = t(:);
    y = y(:);
    p = p(:)';
    if nargin>3 && ~isempty(m)
        z = [ones(length(t),1) ((p(1)*sqrt(p(1)/p(2)^2)./(p(1)+t)./sqrt(p(1)/p(2)^2+t))*ones(1,length(p)-1)).*[ones(length(t),1) cos(t*p(3:(3+m-1))) exp(-t*p(3+m:end))]];
        %z = [ones(length(t),1) ((((1 + t/p(1)).^(-1) .* (1 + p(2)^2*t/p(1))).^(-1/2)) *ones(1,length(p)-1)).*[ones(length(t),1) cos(t*p(3:(3+m-1))) exp(-t*p(3+m:end))]];
    else
        z = [ones(length(t),1) ((p(1)*sqrt(p(1)/p(2)^2)./(p(1)+t)./sqrt(p(1)/p(2)^2+t))*ones(1,length(p)-1)).*[ones(length(t),1) exp(-t*p(3:end))]];
        %z = [ones(length(t),1) ((((1 + t/p(1)).^(-1) .* (1 + p(2)^2*t/p(1))).^(-1/2)) *ones(1,length(p)-1)).*[ones(length(t),1) exp(-t*p(3:end))]];
    end
    if isreal(z)
        [c,resnorm] = lsqnonneg(z,y);
        % c = z\y;
        z = z*c;
        if nargin>4 && ~isempty(bld)
            semilogx(t,y,'o',t,z, 'markersize', 2.5); drawnow;
        end
        %err = sum((y-z).^2./abs(z));
        %err = sum((y-z).^2);
        err = resnorm;
    else
        err = inf;
    end
    
end
