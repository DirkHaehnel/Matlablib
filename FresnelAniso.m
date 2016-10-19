function [rp,rs,tp,ts] = FresnelAniso(q,n1,n2);

% q=n1(1)*n1(2)*tan(theta)./sqrt(n1(2)^2+n1(1)^2*tan(theta).^2);

n1o = n1(1);
n1e = n1(2);
n2o = n2(1);
n2e = n2(2);
w1o = sqrt(n1o^2-q.^2);
w1e = n1o*sqrt(1-q.^2/n1e^2);    
w2o = sqrt(n2o^2-q.^2);
w2e = n2o*sqrt(1-q.^2/n2e^2);    

rp = (w1e*n2o^2-w2e*n1o^2)./(w1e*n2o^2+w2e*n1o^2);
rs = (w1o-w2o)./(w1o+w2o);
tp = sqrt((1+(n2o^2/n2e^2-1).*(q/n2e).^2)./(1+(n1o^2/n1e^2-1).*(q/n1e).^2)).*2*n1o*n2o.*w1e./(w1e*n2o^2+w2e*n1o^2);
ts = 2*w1o./(w1o+w2o);
