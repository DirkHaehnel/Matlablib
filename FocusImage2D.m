function [fx ex] = FocusImage2D(rho,z,volx,phi,flag)

if nargin<4 || isempty(phi)
    phi = 0;
end

if ndims(phi)>2
    if ndims(volx)>2
        fx = volx(:,:,1);
        ex = volx(:,:,2);
    else
        fx = volx;
        ex = volx;
    end
else
    fx = volx(:,:,1,:);
    for j=2:(size(volx,3)+1)/2
        fx=fx+cos((j-1)*phi)*volx(:,:,j,:)+sin((j-1)*phi)*volx(:,:,(end-1)/2+j,:);
    end
    ex = volx(:,:,1,:);
    for j=2:(size(volx,3)+1)/2
        ex=ex+cos((j-1)*(phi+pi))*volx(:,:,j,:)+sin((j-1)*(phi+pi))*volx(:,:,(end-1)/2+j,:);
    end
    if isreal(volx)
        fx = abs(sum(fx,4));
        ex = abs(sum(ex,4));
    else
        fx = sum(abs(fx).^2,4);
        ex = sum(abs(ex).^2,4);
    end
end

if nargout==0
    if nargin>4 && ~isempty(flag)
        if strcmp(flag,'log')
            mpcolor([-flipud(rho);rho],[z;z],log10([flipud(ex);fx]))
            set(gca,'xticklabel',abs(get(gca,'xtick')))
        end
        if strcmp(flag,'surf')
            surf([-flipud(rho);rho],[z;z],[flipud(ex);fx])
            set(gca,'xticklabel',abs(get(gca,'xtick')))
        end
        if strcmp(flag,'horizontal')
            mpcolor([z;z],[-flipud(rho);rho],[flipud(ex);fx])
            set(gca,'yticklabel',abs(get(gca,'ytick')))
        end
        if strcmp(flag,'vertical')
            z0 = min(z(:))-mean(diff(z(1,:)));
            mpcolor([[-flipud(rho);rho] [-flipud(rho);rho]],[2*z0-fliplr([z;z]) [z;z]],[fliplr([flipud(ex);fx]) [flipud(ex);fx]])
            set(gca,'xticklabel',abs(get(gca,'xtick')))
        end
    else
        mpcolor([-flipud(rho);rho],[z;z],[flipud(ex);fx])
        set(gca,'xticklabel',abs(get(gca,'xtick')))
    end
    clear ex fx
end

