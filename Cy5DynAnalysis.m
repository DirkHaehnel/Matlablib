% load D:\Doc\FCS\Cy5Streptavidin\cydynresults.mat
load D:\Joerg\Doc\Gregor\Cy5Dynamics\Cy5dyn.mat dyn
cydyn = dyn; clear dyn;

t = cydyn(2).autotime;
ind = t>=1e-7;
t = t(ind);
clear x
for j=1:length(cydyn)
    x(:,:,:,j) = cydyn(j).auto(ind,1:2,1:2);
end

% for j=6:10
%     x(:,:,:,j-5) = cydyn(j).auto(ind,1:2,1:2);
% end

clear t1 t2
for j=21:size(x,4)
    err = 0;
    ind=t>1e-6 & t<1e-4;
    t2(:,j) = simplex('expfun',[1e-4 1e-3],[],[],[],[],t(ind),1-x(ind,2,1,j),1);
    t1(:,j) = simplex('expfun',[1e-4 1e-3],[],[],[],[],t(ind),1-x(ind,1,2,j),1);
end

clear p c
for j=21:size(x,4)
    err = 0;
    %     ind=t<1e-4;
    %     t2 = simplex('expfun',1e-5,[],[],[],[],t(ind),1-x(ind,2,1,j),1);
    % 	t1 = simplex('expfun',1e-5,[],[],[],[],t(ind),1-x(ind,1,2,j),1);
    ind=t>=1e-7 & t<=1e-3; % | t>=1;
    r = simplex('rigler',[1e-2 1e-3 1e6 1e5],[],[],[],[],t(ind),x(ind,1,1,j)+x(ind,2,1,j)+x(ind,1,2,j)+x(ind,2,2,j))
    %tmp = [r(1:2)' 1./mean(r(3:4)) t1 t2]';
    tmp = [r(1:2)' 1e-5 5e-4 5e-4]';
    for casc=1:10
        p0 = tmp(:,err==min(err));
        for sub=1:10
            tmp(:,sub) = simplex('FLCS',p0.*exp(log(2)*(2*rand(size(p0))-1)),zeros(size(p0)),[],[],[],t(ind),x(ind,:,:,j));
            err(sub) = FLCS(tmp(:,sub),t(ind),x(ind,:,:,j),1);
            drawnow;
        end
    end
    p(:,j) = tmp(:,err==min(err));
    [err, c(:,:,j)] = FLCS(p(:,j),t(ind),x(ind,:,:,j));
end

for j=1:size(x,4)
    ind=t>=1e-6;
    r = simplex('rigler',[1e-2 1e-3 1e6 1e4],[],[],[],[],t(ind),x(ind,1,1,j)+x(ind,2,1,j)+x(ind,1,2,j)+x(ind,2,2,j),[],1)
    p12(:,j) = Simplex('Rigler',[r(1:2)' 1e-4 1e-4],[r(1:2)' zeros(1,2)],[r(1:2)' inf inf],[],[],t(ind),x(ind,1,2,j),[],1);
    [err, c12(:,:,j)] = Rigler(p(:,j),t(ind),x(ind,1,2,j));
end


% save Cy5DynAnalsyis p c name

plot(pow([1 2 5]),1./p(5,[1 2 5]),'o',pow([1 2 5]),polyval(polyfit(pow([1 2 5]),1./p(5,[1 2 5]),1),pow([1 2 5])))
xlabel('power [\muW]')
ylabel('transition rate \itk_{sn}\rm + \itk_{ns}\rm [1/s]')

for j=1:5 r(j) = sqrt(c(5,2,j)/c(5,3,j)/c(5,4,j)*c(5,1,j)); end
plot(pow([1 2 5]),r([1 2 5]),'o',0:400,ones(size(0:400))*mean(r([1 2 5])))
axis([0 400 0 4])
xlabel('power [\muW]')
ylabel('transition rate ratio \itk_{sn}\rm/\itk_{ns}\rm')

ksn = r([1 2 5])./(1+r([1 2 5]))./p(5,[1 2 5]);
kns = 1./(1+r([1 2 5]))./p(5,[1 2 5]);
ksn400 = polyval(polyfit(pow([1 2 5]),ksn,1),0:400);
kns400 = polyval(polyfit(pow([1 2 5]),kns,1),0:400);
ksn400a = polyval(polyfit([-pow([1 2 5]) pow([1 2 5])],[-ksn ksn],1),0:400);
kns400a = polyval(polyfit([-pow([1 2 5]) pow([1 2 5])],[-kns kns],1),0:400);

plot(0:400,ksn400,'r',0:400,kns400,'b',pow([1 2 5]),ksn,'or',pow([1 2 5]),kns,'ob')
xlabel('power [\muW]')
ylabel('transition rate [1/s]')
legend({'\itk_{sn}','\itk_{ns}'},2)

plot(0:400,ksn400,'r',0:400,kns400,'b',0:400,ksn400a,'r:',0:400,kns400a,'b:',pow([1 2 5]),ksn,'or',pow([1 2 5]),kns,'ob')
legend({'\itk_{sn}','\itk_{ns}'},2)
xlabel('power [\muW]')
ylabel('transition rate [1/s]')
