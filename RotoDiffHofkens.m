function [coef, om] = RotoDiffHofkens(drot,nL)

% nL is the max number of angular components
nn = 10; % max size of matrix 

cL0 = (2*(0:nL)+1)/2;

dm = drot(1)-drot(2);
dp = drot(1)+drot(2);
dc = drot(3);

for L=0:nL
    A = zeros(nn); 
    for j=1:nn
        A(j,j) = dp*L*(L+1)/2 + dc*(2*j-2)^2;
        if j<nn
            A(j,j+1) = -dm*L*(L+1)/4;
        end
        if j>1 
            A(j,j-1) = -dm*L*(L+1)/4;
        else
            A(j,j+1) = -dm*L*(L+1)/2;
        end
    end
    [v,d] = eig(A);
    vi = inv(v);
    coef(L+1,:) = vi(:,1)'.*v(1,:)*cL0(L+1);
    om(L+1,:) = diag(d)';
end

return

clear all

nL = 50;
theta = (0.5:1e3)/1e3*pi;
t=0.025*exp(log(2/0.025)/20*(0:20));
for j=0:nL x(j+1,:) = polyval(LegendrePoly(j),cos(theta)); end

drot = [1 1 1];
[coef, om] = RotoDiffHofkens(drot,nL);
for j=1:length(t) z0(j,:)=sum((sum(coef.*exp(-t(j)*om),2)*ones(1,length(theta))).*x); end

drot = [0.5 1.5 1];
[coef, om] = RotoDiffHofkens(drot,nL);
for j=1:length(t) z(j,:)=sum((sum(coef.*exp(-t(j)*om),2)*ones(1,length(theta))).*x); end
plot(theta,z0,theta,z)