function [err,c,z] = RotoDiffFit(p,kv,mm)

kv = kv(:); p = p(:)';

if length(p)==2
    z(:,1) = exp(-(p(1)+p(2))*kv);
    if size(mm,2)>1
        z(:,2) = 1/3 + exp(-6*p(1)*kv)/6 + exp(-2*(p(1)+2*p(2))*kv)/2;
    end
    if size(mm,2)>2
        z(:,3) = 3*exp(-(p(1)+p(2))*kv)/5 + 3*exp(-(11*p(1)+p(2))*kv)/20 + exp(-(3*p(1)+9*p(2))*kv)/4;
    end
    if size(mm,2)>3
        z(:,4) = 1/5 + exp(-6*p(1)*kv)/7 + 9*exp(-20*p(1)*kv)/280 + 3*exp(-(2*p(1)+4*p(2))*kv)/7 + exp(-(16*p(1)+4*p(2))*kv)/14  + exp(-(4*p(1)+16*p(2))*kv)/8;
    end
else
    for k=1:length(p)/2
        z(:,1,k) = exp(-(p(2*k-1)+p(2*k))*kv);
        if size(mm,2)>1
            z(:,2,k) = 1/3 + exp(-6*p(2*k-1)*kv)/6 + exp(-2*(p(2*k-1)+2*p(2*k))*kv)/2;
        end
        if size(mm,2)>2
            z(:,3,k) = 3*exp(-(p(2*k-1)+p(2*k))*kv)/5 + 3*exp(-(11*p(2*k-1)+p(2*k))*kv)/20 + exp(-(3*p(2*k-1)+9*p(2*k))*kv)/4;
        end
        if size(mm,2)>3
            z(:,4,k) = 1/5 + exp(-6*p(2*k-1)*kv)/7 + 9*exp(-20*p(2*k-1)*kv)/280 + 3*exp(-(2*p(2*k-1)+4*p(2*k))*kv)/7 + ...
                exp(-(16*p(2*k-1)+4*p(2*k))*kv)/14  + exp(-(4*p(2*k-1)+16*p(2*k))*kv)/8;
        end
    end
end

if length(p)>2
    col = z(:,:,1); col = col(:);
    for j=2:length(p)/2
        tmp = z(:,:,j);
        col = [col tmp(:)];
    end
    c = col\mm(:);
    tmp = c(1)*z(:,:,1);
    for j=2:length(p)/2
        tmp = tmp + c(j)*z(:,:,j);
    end
    z = tmp;
else
    %     for j=1:size(z,2)
    %         z(:,j) = (z(:,j)\mm(:,j))*z(:,j);
    %     end
    z = (z(:)\mm(:))*z;
end

plot(kv,mm,'o',kv,z); drawnow
err = sum(sum((mm-z).^2./abs(z)));
