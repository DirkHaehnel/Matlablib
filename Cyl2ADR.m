function [adr, em] = Cyl2ADR(theta,phi,coef,mv,n)

theta = theta(:); phi = phi(:)';
st = sin(theta);
em = zeros(length(theta),length(phi),2);
for j=1:length(mv)
    em(:,:,1) = em(:,:,1) + (st.*coef(:,j,1))*exp(i*mv(j)*(phi-pi/2));
    em(:,:,2) = em(:,:,2) + n*(st.*coef(:,j,2))*exp(i*mv(j)*(phi-pi/2));    
end
em = sqrt(n/2/pi)*n*em;
adr = squeeze((abs(em(:,:,1)).^2 + abs(em(:,:,2)).^2));
