function err = BeadScanCorrelationFun(p,t,y);

t = t(:);
z = exp(-t.^2/2/p(1).^2);
z = [z exp(-(t-p(2)).^2/2/p(1).^2) exp(-(t+p(2)).^2/2/p(1).^2)];
col = ones(size(t));

c1 = [col z(:,1)]\y(:,1); 
c2 = [col z(:,1)]\y(:,2); 
% c3 = [col\(y(:,3)-z(:,2)*sqrt(c1(2).*c2(2))); sqrt(c1(2).*c2(2))];
% c4 = [col\(y(:,4)-z(:,3)*sqrt(c1(2).*c2(2))); sqrt(c1(2).*c2(2))];
c3 = [col z(:,2)]\y(:,3); 
c4 = [col z(:,3)]\y(:,4); 

z = [[col z(:,1)]*c1 [col z(:,1)]*c2 [col z(:,2)]*c3 [col z(:,3)]*c4];

plot(t,y,'o',t,z); drawnow

err = sum(sum((y-z).^2));


