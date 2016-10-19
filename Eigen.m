% program Eigen

if true
    
    bj = inline('besselj(n+0.5,x).*sqrt(pi./x/2)','n','x');
    by = inline('bessely(n+0.5,x).*sqrt(pi./x/2)','n','x');
    bh = inline('besselh(n+0.5,x).*sqrt(pi./x/2)','n','x');
    bjprime = inline('(besselj(n-0.5,x)-besselj(n+1.5,x)-besselj(n+0.5,x)./x).*sqrt(pi./x/2)/2','n','x');
    byprime = inline('(bessely(n-0.5,x)-bessely(n+1.5,x)-bessely(n+0.5,x)./x).*sqrt(pi./x/2)/2','n','x');
    bhprime = inline('(besselh(n-0.5,x)-besselh(n+1.5,x)-besselh(n+0.5,x)./x).*sqrt(pi./x/2)/2','n','x');
    
    clear z zz
    dq = 0.01; q = (dq/2:dq:100)';
    % q = 5;
    x1 = 1.025:0.05:2;
    x2 = 2.025:0.05:3;
    x = [x1 x2];
    row1 = ones(1,length(x1));
    row2 = ones(1,length(x2));
    qq1 = q*row1;
    qq2 = q*row2;
    
    
    %     load bjzeros
    %     t = 0.1; 
    %     z = zeros(length(x),41);
    %     for k=0:50
    %         q = rb(k+1,:)';
    %         z(:,k+1) = sum(bj(k,q*x).*((q.^2.*by(k,q*2).*exp(-t*q.^2))*ones(1,length(x))))'; 
    %         k
    %     end
    %     zz = zeros(size(z,1),1); for j=1:size(z,2) zz=zz+(2*j-1)*z(:,j); end
    %     plot(x,zz); figure(gcf)
    
    
    t = 0.1; 
    z = zeros(length(x),41); y = z;
    for k=0:75
        z(:,k+1) = sum(bj(k,q*x).*((q.^2.*bj(k,q*2).*exp(-t*q.^2))*ones(1,length(x))))'; 
        y(:,k+1) = sum((by(k,q*x)./((by(k,q)./bj(k,q))*ones(size(x)))).*((q.^2.*bj(k,q*2).*exp(-t*q.^2))*ones(1,length(x))))'; 
        k
    end
    zz = zeros(size(z,1),1); for j=1:size(z,2) zz=zz+(2*j-1)*z(:,j); end
    zz = zz*dq;
    yy = zeros(size(y,1),1); for j=1:size(y,2) yy=yy+(2*j-1)*y(:,j); end
    yy = yy*dq;
    plot(x,zz,'.',x,exp(-(x-2).^2/4/t)/sqrt(4/pi*t^3)/2); figure(gcf)
    
    
    % t = 0.02; 
    % z = zeros(length(x),41);
    % for k=0:40
    %     z(1:length(x1),k+1) = sum(bj(k,q*x1).*((q.^2.*bj(k,q*2))*row1).*exp(-t*qq1.^2))'; 
    %     z(length(x1)+(1:length(x2)),k+1) = sum(bj(k,q*x2).*((q.^2.*bj(k,q*2))*row2).*exp(-t*qq2.^2))';
    %     k
    % end
    % zz = zeros(size(z,1),1); for j=1:size(z,2) zz=zz+(2*j-1)*z(:,j); end
    % zz = zz*dq;
    % plot(x,zz,'.',x,exp(-(x-2).^2/4/t)/sqrt(4/pi*t^3)/2); figure(gcf)
    
    
    % t = 0.1; 
    % for k=0:100
    %     tmp1 = bj(k,q*x1); tmp2 = bh(k,q*x2);
    %     z(1:length(x1),k+1) = (sum(((tmp1 + (q*x1).*bjprime(k,q*x1)).*(bh(k,q*2)*row1)+2*(q*row1).*tmp1.*bhprime(k,qq1*2)).*exp(-t*qq1.^2)))'; 
    %     z(length(x1)+(1:length(x2)),k+1) = ...
    %         (sum(((tmp2 + (q*x2).*bhprime(k,q*x2)).*(bj(k,q*2)*row2)+2*(q*row2).*tmp2.*bjprime(k,qq2*2)).*exp(-t*qq2.^2)))';
    %     k
    % end
    % zz = zeros(size(z,1),1); for j=1:size(z,2) zz=zz+(2*j-1)*z(:,j); end
    % zz = zz*dq;
    % plot(x,zz,'.',x,exp(-(x-2).^2/4/t)/sqrt(4/pi*t)); figure(gcf)
    
    % for k=0:100
    %     z(1:length(x1),k+1) = bj(k,q*x1')*(bh(k,q*2) + q*2*bhprime(k,q*2)) + q*x1'.*bjprime(k,q*x1')*bh(k,q*2); 
    %     z(length(x1)+(1:length(x2)),k+1) = (bj(k,q*2) + q*2*bjprime(k,q*2)).*bh(k,q*x2') + q*bj(k,q*2)*x2'.*bhprime(k,q*x2');
    %     k
    % end
    % zz = zeros(size(z,1),1); for j=1:size(z,2) zz=zz+(2*j-1)*z(:,j); end
    % plot(x,zz,'.',x,real(exp(i*q*(x-2)))); figure(gcf)
    
    % for k=0:100
    %     z(1:length(x1),k+1) = bj(k,q*x1').*bh(k,q*2); 
    %     z(length(x1)+(1:length(x2)),k+1) = bj(k,q*2).*bh(k,q*x2');
    %     k
    % end
    % zz = zeros(size(z,1),1); for j=1:size(z,2) zz=zz+(2*j-1)*z(:,j); end
    % zz = i*q*zz; plot(x,zz,'.',x,real(exp(i*q*(x-2))./abs(x-2))); figure(gcf)
    
end


if false

    dist = 50e-7; % distance between receptor and source
    kp = 1e8; % on rate in 1/M^2/s
    km = 1e2; % off rate in 1/M/s
    diffconst = 1e-6; % cm^2/s
    dt = 1e-7;
    maxt = 100e-6; 
    t = (dt/2:dt:maxt);
    r = (0.5:100)'/100*dist;
    col = ones(size(r));
    
    c = exp(-r.^2*(1./(4*diffconst*t)))./(4*pi*diffconst*col*t).^(3/2);
end