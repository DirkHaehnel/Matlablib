load D:\Joerg\Doc\IngoW\cGMP_IBMX.mat 
% load D:\Doc\IngoW\cGMP_IBMX.mat 

s = int2str(4);

eval(['sperm = spermv( ' s ')*1e3;']);
eval(['t = t' s]);
eval(['conc = union(c' s ',c' s ');']);
clear y
for j=1:length(conc)
    eval(['tmp = x' s '(:,c' s '==conc(j));']);
    ind = isnan(tmp);
    tmp(ind) = 0;
    y(:,j) = sum(tmp,2)./sum(double(~ind),2);
end
for j=1:size(y,2) y(:,j) = abs(y(:,j)-y(1,j)); end
%ind = 1:length(t)-1;
ind = 1:length(t);

tmp = [1.6e-7 2.8e7 1.13e6 270 0.5]';
for casc=1:10
    tmp = tmp(:,err==min(err));
    err = KinFun(tmp,conc,t(ind),y(ind,:),1);
    for sub=1:3
        p = tmp(:,1).*exp(log(2)*(2*rand(size(p0))-1))
        tmp(:,sub+1) = simplex('KinFun',p,0*p,[],[],[],conc,t(ind),y(ind,:),1);
        err(sub+1) = KinFun(tmp(:,sub+1),conc,t(ind),y(ind,:),1);
        drawnow;
    end
end
p = tmp(:,err==min(err));

[err,c] = KinFun(p,conc,t(ind),y(ind,:),1);
s = [num2str(1e9*conc',3) repmat(' nM',length(conc),1)];
h = legend(s,2); set(h,'fontsize',12)
xlabel('time [s]')
ylabel('cGMP [pmol / 10^8 sperm]');
pp = p;
pp(1) = pp(1)*6.0221e23/sperm;
pp(3) = c*sperm/1e8*1e-12*pp(3);
pp(4) = pp(4)*(sperm/1e8*1e-12/c)^pp(5);
gtext({['receptors/sperm = ' num2str(pp(1),3)],...
    ['binding constant = ' num2str(pp(2),3) ' M^{-1}\cdots^{-1}'],...
    ['guanylate cyclase activity = ' num2str(pp(3),3) ' s^{-1}'],...
    ['guanylate cyclase inhibition = ' num2str(pp(4),3) ' M^{' num2str(-pp(5),3) '}\cdots^{-1}'],...
    ['inhibition coefficient = ' num2str(pp(5),3)]},'fontsize',12)

print -dpng res6
save res6 p err


return

p = simplex('KinFun',p,0*p,[],[],[],conc,t(ind),y(ind,:),1);
err = KinFun(p,conc,t(ind),y(ind,:),1);
for sub=1:500
    p(:,end+1) = simplex('KinFun',p(:,end),0*p(:,1),[],[],[],conc,t(ind),y(ind,:),1);
    err(end+1) = KinFun(p(:,end),conc,t(ind),y(ind,:),1);
    drawnow;
end
p(:,err==min(err))


return

err = [];
vv = [10.^-(3:6); 10.^(4:2:10); 10.^(-2:2:4); 10.^(2:2:8); 10.^(3:6); 0.6:0.2:1.2]';
for j1=1:4
    for j2=1:4
        for j3=1:4
            for j4=1:4
                for j5=1:4
                    for j6=1:4
                        err(j1,j2,j3,j4,j5,j6) = KinFun([vv(j1,1) vv(j2,2) vv(j3,3) vv(j4,4) vv(j5,5) vv(j6,6)],...
                            conc,t(1:end-1),y(1:end-1,:));
                    end
                end
            end
        end
    end
end
[k1,k2,k3,k4,k5,k6]=ndgrid(1:4,1:4,1:4,1:4,1:4,1:4);
[k1(err==min(err(:))) k2(err==min(err(:))) k3(err==min(err(:))) k4(err==min(err(:))) k5(err==min(err(:))) k6(err==min(err(:)))]

tmp = [vv(k1(err==min(err(:))),1) vv(k2(err==min(err(:))),2) vv(k3(err==min(err(:))),3) vv(k4(err==min(err(:))),4) vv(k5(err==min(err(:))),5) vv(k6(err==min(err(:))),6)]'
