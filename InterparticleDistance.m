% numerical simulation of the average nearest neighbor distance of randomly distributed particles in 2D

res = [];
len = 3e3;
for j=1:100
    x = sqrt(len)*rand(1,len);
    y = sqrt(len)*rand(1,len);
    r = sqrt((ones(len,1)*x-x'*ones(1,len)).^2+(ones(len,1)*y-y'*ones(1,len)).^2);
    r(r==0) = inf;
    res = [res mean(min(r))];
end
mean(res)