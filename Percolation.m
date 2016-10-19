close all
rmax = 500;
res = zeros(9,1);
for k=1:9 
   for j=1:100
      [t, len] = cluster(round((0.5+k*0.1)*rand(100)));
      tmp = hist(len,1:max(len))
      if max(len)>size(res,2)
         res = [res zeros(size(res,1), max(len)-size(res,2))];
      end
      res(k,1:max(len)) = res(k,1:max(len)) + tmp;
      semilogy(res(k,:)); drawnow
   end
end

