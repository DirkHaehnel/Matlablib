function p = AntibunchRotoDiffFit(autotime,auto,p0)

t = [-autotime(end:-1:1), autotime];
y = [[auto(end:-1:1,1,4,1);auto(:,4,1,1)],...
    [auto(end:-1:1,2,3,1);auto(:,3,2,1)]];

period = 12.5;
dt = mean(diff(autotime));
n = round(period/dt);

tst = conv(y(:,1),y(:,1));
[ind,ind] = max(tst);
ind = round((ind - size(y,1))/2);
t = t(1:end-ind+1);
y = y(ind:end,:);



p = Simplex('ExpFun', p0, zeros(size(p0)), [], [], [], t, z, 1);

return

p(:,2) = Simplex('ExpFun', p0, [], [], [], [], [-t3(end:-1:1),t4], [z3(end:-1:1) z4], 1);

return


t = [-autotime(end:-1:1), autotime];
y = [[sum(auto(end:-1:1,1,4,:),4);sum(auto(:,4,1,:),4)],...
    [sum(auto(end:-1:1,2,3,:),4);sum(auto(:,3,2,:),4)]];

period = 12.5;
dt = mean(diff(autotime));
n = round(period/dt);

tst = 0;
for j=1:1.5*n
    len = floor((size(y,1)-j)/n)*n;
    tmp = sum(reshape(y(j:j+len-1,1),n,len/n));
    %plot(tmp); drawnow
    tmp = max(tmp)-min(tmp);
    if tmp>tst
        jj = j;
        tst = tmp;
    end
end
j = jj;
len = floor((size(y,1)-j)/n)*n;
t = mean(reshape(t(j:j+len-1),n,len/n));
z(1,:) = sum(reshape(y(j:j+len-1,1),n,len/n));
z(2,:) = sum(reshape(y(j:j+len-1,2),n,len/n));
[ind, ind] = min(z(1,:));
t = t-t(ind);
plot(t,z)

p = Simplex('ExpFun', p0, [], [], [], [], t, z(1,:), 1);